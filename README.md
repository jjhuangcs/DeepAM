# [Learning Deep Analysis Dictionaries for Image Super-Resolution](https://ieeexplore.ieee.org/document/9257106)

<!--![visitors](https://visitor-badge.glitch.me/badge?page_id=jjhuangcs/WINNet)-->

[Jun-Jie Huang](https://jjhuangcs.github.io/) (jjhuang@nudt.edu.cn) and [Pier Luigi Dragotti](http://www.commsp.ee.ic.ac.uk/~pld/)

Matlab implementation for "Learning Deep Analysis Dictionaries for Image Super-Resolution" (TSP'2020).

<img width="654" alt="DeepAM" src="https://user-images.githubusercontent.com/89965355/196571760-f3f7f6a4-5473-46f7-bf0e-e0e0dcf65e1f.png">


Inspired by the recent success of deep neural networks and the recent efforts to develop multi-layer dictionary models, we propose a Deep Analysis dictionary Model (DeepAM) which is optimized to address a specific regression task known as single image super-resolution. Contrary to other multi-layer dictionary models, our architecture contains L layers of analysis dictionary and soft-thresholding operators to gradually extract high-level features and a layer of synthesis dictionary which is designed to optimize the regression task at hand. In our approach, each analysis dictionary is partitioned into two sub-dictionaries: an Information Preserving Analysis Dictionary (IPAD) and a Clustering Analysis Dictionary (CAD). The IPAD together with the corresponding soft-thresholds is designed to pass the key information from the previous layer to the next layer, while the CAD together with the corresponding soft-thresholding operator is designed to produce a sparse feature representation of its input data that facilitates discrimination of key features. DeepAM uses both supervised and unsupervised setup. Simulation results show that the proposed deep analysis dictionary model achieves better performance compared to a deep neural network that has the same structure and is optimized using back-propagation when training datasets are small. On noisy image super-resolution, DeepAM can be well adapted to unseen testing noise levels by rescaling the IPAD and CAD thresholds of the first layer.

# How to use the code?
1. use trainDeepAM.m to perform training.

2. use testDeepAM.m to perform testing on Set5 or Set14. 

In the folder models, there are 3 learned models using trainDeepAM.m corresponding to learning without noise and with noise level 0.1.


# Citation

If you use any part of this code in your research, please cite our paper:


```
@ARTICLE{Huang2020DeepAM,
	author={Huang, Jun-Jie and Dragotti, Pier Luigi},
  	journal={IEEE Transactions on Signal Processing}, 
  	title={Learning Deep Analysis Dictionaries for Image Super-Resolution}, 
  	year={2020},
  	volume={68},
  	number={},
  	pages={6633-6648},
  	doi={10.1109/TSP.2020.3036902}
}
```

