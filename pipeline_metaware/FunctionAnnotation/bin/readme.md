# 宏基因组功能注释模块

## KEGG 各个id间功能关系数据库构建 

- api list: 获取全部的KO
- api link: 基于KO,获取对应的pathway（一个KO可能对应多个pathway）
- api get: 基于map-id, 获取level信息

### 注意
下载的level文件中,有部分不是以map开头的, 具体原因详见：

```
 the following meaning:

map
    manually drawn reference pathway
ko
    reference pathway highlighting KOs
ec
    reference metabolic pathway highlighting EC numbers
rn
    reference metabolic pathway highlighting reactions
<org>
    organism-specific pathway generated by converting KOs to gene identifier

```