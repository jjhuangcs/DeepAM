# [Learning Deep Analysis Dictionaries for Image Super-Resolution](https://ieeexplore.ieee.org/document/9257106)

<!--![visitors](https://visitor-badge.glitch.me/badge?page_id=jjhuangcs/WINNet)-->

[Jun-Jie Huang](https://jjhuangcs.github.io/) (jjhuang@nudt.edu.cn) and [Pier Luigi Dragotti](http://www.commsp.ee.ic.ac.uk/~pld/)

Matlab implementation for "Learning Deep Analysis Dictionaries for Image Super-Resolution" (TSP'2020).

<img width="654" alt="WINNet" src="https://user-images.githubusercontent.com/89965355/178172283-b6b9e7da-add2-44ad-b83d-3b87918a8c5b.png">

#How to use the code?
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

