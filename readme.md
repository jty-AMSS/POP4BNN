# Verifying Properties of Binary Neural Networks Using Sparse Polynomial Optimization

This repository is the official implementation for "Verifying Properties of Binary Neural Networks Using Sparse Polynomial Optimization". 

## Requirements

For Julia, it is necessary to install the following packages: TSSOS, MosekTools (needs license), DynamicPolynomials and DelimitedFiles;
In Julia, run 
```
pkg> add https://github.com/wangjie212/TSSOS
pkg> add DelimitedFiles
```

For Matlab, it is necessary to install  the Gurobi (needs license). 

For Python (to train BNN),  it is necessary to install   tensorflow  and larq, run 

```
pip install tensorflow 
pip install larq
```


## Evaluation



To verify  Robustness against $\| \cdot \|_{\infty}-$ attacks by our method, run  'Exp_SDP.m' in Matlab.
To verify  Robustness against $\| \cdot \|_{\infty}-$ attacks by Gurobi, run  'Exp_Grb.m' in Matlab.

To verify  Robustness against $\| \cdot \|_{2}-$ attacks by our method, run  'Exp_SDP_L2_with_Bound.m' in Matlab.
To verify  Robustness against $\| \cdot \|_{2}-$ attacks by Gurobi, run  'Exp_Grb_L2_with_bound.m' in Matlab.

Please modify the parameters 'NetDataPath' and 'eps_tol' in the corresponding '.m' file to compute the results for different BNNs and different values of $\varepsilon$.


For the source code of the standard verification problem, please refer to 'BNNV2POP_Large.m'.


## Training

To train the BNNs and get the parameters in JSON format, run this command:

```train
python TrainBNN.py
```


## Results



### [Robustness against $\| \cdot \|_{\infty}-$ attacks]


For $\text{BNN}_1:$ with size $[784, 500, 500, 10], w_s=34.34\%$:
| $\delta_{\|\|\cdot\|\|_{\infty}}$ | $ \tau_{\text{LP}}$  verified   |  $ \tau_{\text{LP}}$  time (s)    | $\tau_{\text{tighter}, cs}^{1}$  verified   | $\tau_{\text{tighter}, cs}^{1}$  time (s)  |  $\tau_{\text{Soft-MILP}}$ verified        |  $\tau_{\text{Soft-MILP}}$ time  (s)    |
|  ----  | ----  | ----  | ----  |----  | ----  | ----  |
| $0.25$|       83|     0.01| 91|  3.62  | 3 NR+95 RO| 0.04|
| $0.50$|      31|    0.02| 60|  6.69  | 4 NR+94 RO| 1.21|
| $1.00$|     1|    0.03| 21|  10.76|  15 NR+50 RO|  251.90|%+33 TO 
| $1.50$|       0|   0.06|  6|  38.32|  20 NR+12 RO| 428.24|  


For $\text{BNN}_2:$ with size $[784, 500, 500, 10], w_s=19.07\%$:
| $\delta_{\|\|\cdot\|\|_{\infty}}$ | $ \tau_{\text{LP}}$  verified   |  $ \tau_{\text{LP}}$  time (s)    | $\tau_{\text{tighter}, cs}^{1}$  verified   | $\tau_{\text{tighter}, cs}^{1}$  time (s)  |  $\tau_{\text{Soft-MILP}}$ verified        |  $\tau_{\text{Soft-MILP}}$ time  (s)    |
|  ----  | ----  | ----  | ----  |----  | ----  | ----  |
|$0.25$|     14 |     0.03 | 59 | 11.97  | 3 NR+95 RO|2.23|
|$0.50$|      0 |    0.05 | 23 |  42.37|  9 NR+63 RO| 220.53|
|$0.75$|      0 |    0.08 | 9 |  139.18|  10 NR+19 RO|  455.61|
### [Robustness against $\| \cdot \|_{2}-$ attacks]

For $\text{BNN}_1:$ with size $[784, 500, 500, 10], w_s=34.34\%$:


| $\delta_{\|\|\cdot\|\|_{2}}$ | $\tau_{\text{tighter}, cs}^{1}$   verified  | $\tau_{\text{tighter}, cs}^{1}$  time (s)  |  $\tau_{\text{Soft-MINP}}$ verified        |  $\tau_{\text{Soft-MINP}}$ time  (s)    |
|  ----  | ----  | ----  | ----  |----  |  
|$10$|      70 |5.23 |  3 NR+93 RO |  33.35 |
|$20$|     36  |  19.54  |4 NR+30 RO |   447.11 |
|$30$|     13 |   34.24 |  4 NR+ 6  RO |   556.07|

For $\text{BNN}_2:$ with size $[784, 500, 500, 10], w_s=19.07\%$:
| $\delta_{\|\|\cdot\|\|_{2}}$ | $\tau_{\text{tighter}, cs}^{1}$   verified  | $\tau_{\text{tighter}, cs}^{1}$  time (s)  |  $\tau_{\text{Soft-MINP}}$ verified        |  $\tau_{\text{Soft-MINP}}$ time  (s)    |
|  ----  | ----  | ----  | ----  |----  |  
|$5$|     81 |8.57 | 2 NR+96 RO |  3.31 |
|$10$|     46 |19.00 | 3 NR+ 58 RO  |  272.92 |
|$15$|     27 |63.21 | 4 NR+ 23 RO |  475.78 |
