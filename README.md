# Dimensionality Reduction for Regression

An unsupervised method for dimensionality reduction. To be implemented in Python using the scikit-learn [framework](http://scikit-learn.org/stable/developers/contributing.html#rolling-your-own-estimator).


### Resources

* [Original Webpage](http://isp.uv.es/drr.html)
* [Original Code](http://isp.uv.es/code/featureextraction/DRR_toolbox_v1.zip)
* [Paper](https://arxiv.org/abs/1602.00214) - Dimensionality Reduction via Regression in Hyperspectral Imagery

#### Abstract from Paper

> This paper introduces a new unsupervised method for dimensionality reduction via regression (DRR). The algorithm belongs to the family of invertible transforms that generalize Principal Component Analysis (PCA) by using curvilinear instead of linear features. DRR identifies the nonlinear features through multivariate regression to ensure the reduction in redundancy between he PCA coefficients, the reduction of the variance of the scores, and the reduction in the reconstruction error. More importantly, unlike other nonlinear dimensionality reduction methods, the invertibility, volume-preservation, and straightforward out-of-sample extension, makes DRR interpretable and easy to apply. The properties of DRR enable learning a more broader class of data manifolds than the recently proposed Non-linear Principal Components Analysis (NLPCA) and Principal Polynomial Analysis (PPA). We illustrate the performance of the representation in reducing the dimensionality of remote sensing data. In particular, we tackle two common problems: processing very high dimensional spectral information such as in hyperspectral image sounding data, and dealing with spatial-spectral image patches of multispectral images. Both settings pose collinearity and ill-determination problems. Evaluation of the expressive power of the features is assessed in terms of truncation error, estimating atmospheric variables, and surface land cover classification error. Results show that DRR outperforms linear PCA and recently proposed invertible extensions based on neural networks (NLPCA) and univariate regressions (PPA). 
