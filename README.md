# drfDLNMadditive
R package for additive Distributed Response Function Distributed Lag Non-Linear Model (Additive DRF-DLNM), which is proposed in the paper:

> Gasparrini, A., Scheipl, F., Armstrong, B., & Kenward, M. G. (2017). A penalized framework for distributed lag non-linear models. *Biometrics*, 73(3), 938-948.

We call this model "DRF-DLNM" to distinguish it from the Adaptive Cumulative Exposure DLNM (ACE-DLNM) model proposed in: 

> Pan, T., Shin, H. H., McGee, G., & Stringer, A. (2025). Estimating associations between cumulative exposure and health via generalized distributed lag non-linear models using penalized splines. *Biometrics*, 81(3), ujaf116.

Instead of using `mgcv`, we implement the additive DRF-DLNM in `Cpp`. The implementation supports the fast calculation of leave-one-out cross validation. 





## Installation
+ Install the development version from GitHub:

```R
# install.packages('devtools')
devtools::install_github("tianyi-pan/drfDLNMadditive")
```