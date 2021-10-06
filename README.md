# Bioencoder
> 一个面向机器学习的氨基酸序列编码工具

## 简介
一个面向机器学习的氨基酸序列编码工具，主要有如下特点：
- 面向机器学习，返回numpy类型数据，无需处理
- 丰富的编码种类（NUM、BE、EAAC、AAINDEX、GACC、CKSAAP）
- Fasta 原生支持，支持超大文件读取
- 开箱即用（pip 安装方式 简单易用）
  
## 安装方法
```bash
python setup.py sdist
pip install ./dist/bioencoder-beta.tar.gz
```
## 使用方法
```python
from bioencoder import *
pos_data = "pos.fasta"
window_size = 12
pos_seqList,pos_labellist,pos_seqNamelist=get_data(pos_data,1,method="GAAC",window_size=window_size)
```
## 作者
张淋

