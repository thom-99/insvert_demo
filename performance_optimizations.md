The proposed fix targets file I/O by batching fh.write() calls into 1MB blocks, but benchmarking shows this has virtually no effect (16041ms → 15711ms on a 5MB chromosome) because file I/O is not the bottleneck. The real problem is the string shrinking pattern inside the write loop — self.buffer = self.buffer[self.width:] — which on every iteration creates a full copy of the entire remaining string, just 60 characters shorter. For a 5MB chromosome that means ~83,000 full string copies of decreasing size, which is O(n²) memory allocation behaviour. Isolating just the slicing loop with no file I/O at all confirms this: the shrinking pattern takes 50ms on 500KB, while an index-based approach that slices from a fixed string takes 1.4ms — a 35× difference with no I/O involved whatsoever. The fix is to never shrink the string, and instead walk through it with a position index, only copying the small leftover tail at the end:

```
# CURRENT - O(n²): creates a full copy of the remaining string on every iteration
def write(self, data):
    self.buffer += data.upper()
    while len(self.buffer) >= self.width:
        self.fh.write(self.buffer[:self.width] + '\n')
        self.buffer = self.buffer[self.width:]  # ← full string copy each time
```

```
# FIXED - O(n): walks the string with an index, copies only the small leftover
def write(self, data):
    self.buffer += data.upper()
    pos = 0
    lines = []
    while pos + self.width <= len(self.buffer):
        lines.append(self.buffer[pos:pos + self.width])
        pos += self.width
    if lines:
        self.fh.write('\n'.join(lines) + '\n')
    self.buffer = self.buffer[pos:]  # only the tail (<60 chars) is ever copied
```
