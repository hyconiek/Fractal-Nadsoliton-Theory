# FR224 next-atom diagnostic

FR223 best gap: `-8.537388502427579e-06` (negative; navigation only).

Best point:

```text
x=0.00014005601545497015
u=0.00012451172122129872
v=0.00016451613150714195
e=1.0200000366878254e-05
y=1.1276252044436085e-05
```

Navigation buffer: `1.02`.

Nearest buffered shifted-mask walls:

- **FR100**: max relative violation `1.504e-08`; components `{'x': 0, 'u': 0, 'v': 1.504341180428245e-08, 'e': 0}`.
- **FR102**: max relative violation `1.504e-08`; components `{'x': 0, 'u': 0, 'v': 1.504341180428245e-08, 'e': 0}`.
- **FR184**: max relative violation `1.504e-08`; components `{'x': 0, 'u': 0, 'v': 1.504341180428245e-08, 'e': 0}`.
- **FR28**: max relative violation `3.597e-08`; components `{'x': 0, 'u': 0, 'v': 0, 'e': 3.596845623333453e-08}`.
- **FR40**: max relative violation `3.597e-08`; components `{'x': 0, 'u': 0, 'v': 0, 'e': 3.596845623333453e-08}`.
- **FR110**: max relative violation `3.597e-08`; components `{'x': 0, 'u': 0, 'v': 0, 'e': 3.596845623333453e-08}`.
- **FR92**: max relative violation `4.965e-08`; components `{'x': 4.965151324755247e-08, 'u': 0, 'v': 0, 'e': 0}`.
- **FR150**: max relative violation `4.965e-08`; components `{'x': 4.965151324755247e-08, 'u': 0, 'v': 0, 'e': 0}`.
- **FR84**: max relative violation `4.965e-08`; components `{'x': 4.965151324755247e-08, 'u': 0, 'v': 1.504341180428245e-08, 'e': 0}`.
- **FR108**: max relative violation `4.965e-08`; components `{'x': 4.965151324755247e-08, 'u': 0, 'v': 1.504341180428245e-08, 'e': 0}`.
- **FR52**: max relative violation `9.766e-04`; components `{'x': 0, 'u': 0.0009765823673033413, 'v': 0, 'e': 3.596845623333453e-08}`.
- **FR50**: max relative violation `4.839e-02`; components `{'x': 0, 'u': 0, 'v': 0.04838711254551236, 'e': 3.596845623333453e-08}`.

Interpretation:

- The smallest violation is the upper-`v` wall shared by **FR100/FR102/FR184**: the point is only about `1.5e-8` relatively outside the 2%-buffered upper `v` face.
- A nearly comparable wall is the `e` face of FR28/FR40/FR110, but FR184 carries a much stronger parity radius than those old `e<=1e-5` boxes.
- Recommended FR224: first try extending **FR184 upward in `v`** while keeping its `x,u,e` ranges. If a single wide box fails by dependency overestimation, use one adjacent off-centre slab across `v=1/6200`.
- After freezing FR224, rerun the fixed-seed residual search on the exact same compact domain and updated mask union.
- If the residual repeatedly returns to `y=1e-6` rather than a local wall, revisit the FR42 global-tail upgrade. The earlier `delta=1e-5` attempt timed out computationally and was not a mathematical failure.
