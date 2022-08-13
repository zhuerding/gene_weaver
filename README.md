# Gene Weaver
## 基因编织
#### 当前版本：***V1.5.0001***

---
目录
- 一、[项目简介](#index1)
- 二、[功能详解](#index2)  
  [1.GSE数据集导入](#title1)

## <span id='index1'>一、项目简介</span>
  欢迎使用Gene Waver。\
  这是一款基于***python3.9***开发的**完全免费**、**开源**的基因筛查工具，本程序利
  用OCR、计算机爬虫等技术可以实现对StarBase、Ualcan、 mirWalk、mirDIP、miRDB、TargetScan
、TarBase、PubMed、中国知网等多个数据库进行检索，获取相应靶基因信息、差异表达、
生存曲线等信息，从而得出可靠的miRNA -> gene symbol通路。同时程序
还可以生成表格、韦恩图等多种输出形式，亦可在公共服务器上运行。  
  

  **适用人群：** 有大批量miRNA检索需求的生信科研人员
## <span id='index2'>二、功能详解</span>  
- ### <span id='title1'>GSE数据集导入</span>  
本程序提供两种导入GSE数据集方式
1. ini配置文件导入  
  在` ./config.ini`文件中的`[base]`模块中米娜桑可以看到`gse = ['','','']`
这一行属性值，当GSE数据集位于**本程序所在**的文件夹时米娜桑只需要按照`['xxx.tsv','xxxx.tsv']`的格式
添加数据即可。  
  如果米娜桑的数据集不在**本程序所在**的文件夹时就需要填写绝对路径，例如
`['C:\Users\xx\Desktop\xxx.tsv','C:\Users\xx\Desktop\xxx.tsv`

  为避免上一次程序运行时未即时清除配置文件中写入的数据集，故在自动导入时会出现提示，让您确认是否使用该数据集地址。
所以请不要在使用自动导入功能后直接就去摸鱼，不然程序也会摸鱼的。
2. 手动导入  
  当出现`手动导入数据集模式启动：`提示时
就说明配置文件中没有正确导入或者没有写入数据集地址，此时便要手动导入数据集