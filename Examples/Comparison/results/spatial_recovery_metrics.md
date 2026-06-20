# Spatial Deconvolution Recovery

Synthetic benchmark with known cell-type fractions per spot. Lower RMSE/MAE/JSD and higher correlation/dominant accuracy indicate better recovery.

- Genes: 90
- Reference cells: 180
- Spatial spots: 12

| Method | RMSE | MAE | Mean Correlation | Dominant Accuracy | Mean JSD |
|---|---:|---:|---:|---:|---:|
| RCTD-like | 0.0855 | 0.0384 | 0.9737 | 0.9167 | 0.0789 |
| Cell2Location-like | 0.0207 | 0.0140 | 0.9981 | 0.9167 | 0.0474 |

