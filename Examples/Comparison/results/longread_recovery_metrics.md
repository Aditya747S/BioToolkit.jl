# Long-read Structural Variant Recovery

Synthetic benchmark with known DEL/INV/TRA truth events. Metrics are computed by event-type aware matching with +/-70 bp positional tolerance.

- Truth events: 3
- Predicted calls: 4

| SV Type | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| ALL | 3 | 1 | 0 | 0.750 | 1.000 | 0.857 |
| DEL | 1 | 0 | 0 | 1.000 | 1.000 | 1.000 |
| INV | 1 | 1 | 0 | 0.500 | 1.000 | 0.667 |
| TRA | 1 | 0 | 0 | 1.000 | 1.000 | 1.000 |

