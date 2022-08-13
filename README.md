# Gene Weaver
## 基因编织
#### 当前版本：***V1.5.0001***

---
目录
- 一、[项目简介](#index1)
- 二、[文件结构](#index2)
- 三、[功能详解](#index3)  
  [1.GSE数据集导入](#title1)  
  [2.StarBase查询](#title2)  
  [3.Unlcan查询](#title3)  
  [4.中国知网&PubMed查询](#title4)  
  [5.mirWalk查询](#title5)  
  [6.miRDB查询](#title6)  
  [7.TargetScan查询](#title7)  

## <span id='index1'>一、项目简介</span>
  欢迎使用Gene Waver。\
  这是一款基于***python3.9***开发的**完全免费**、**开源**的基因筛查工具，本程序利
  用OCR、计算机爬虫等技术可以实现对StarBase、Ualcan、 mirWalk、mirDIP、miRDB、TargetScan
、TarBase、PubMed、中国知网等多个数据库进行检索，获取相应靶基因信息、差异表达、
生存曲线等信息，从而得出可靠的miRNA -> gene symbol通路。同时程序
还可以生成表格、韦恩图等多种输出形式，亦可在公共服务器上运行。  
  

  **适用人群：** 有大批量miRNA检索需求的生信科研人员  
## <span id='index2'>二、文件结构</span>  
- `./gene_waver.exe`  主程序
- `./database/hsa-lists.csv`  miRNA ID互译表
- `./README.md`  用户手册（markdown语法）
- `./cancer.txt`  癌种代号
- `./test.csv`  测试/体验 GSE数据集
- `./config.ini`  程序配置文件，包含程序运行所需参数  

其余文件为主程序依赖的包，**请勿修改**

## <span id='index3'>三、功能详解</span>  
- ### <span id='title1'>GSE数据集导入</span>  
本程序提供两种导入GSE数据集方式
1. ini配置文件导入  
  在` ./config.ini`文件中的`[base]`模块中米娜桑可以看到`gse = ['','','']`这一行属性值，当GSE数据集位于**本程序所在**的文件夹时米娜桑只需要按照`['xxx.tsv','xxxx.tsv']`的格式 添加数据即可。  
  如果米娜桑的数据集不在**本程序所在**的文件夹时就需要填写绝对路径，例如`['C:\Users\xx\Desktop\xx.tsv','C:\Users\xx\Desktop\xx.tsv`。  
  为避免上一次程序运行时未即时清除配置文件中写入的数据集，故在自动导入时会出现提示，让您确认是否使用该数据集地址。
所以请不要在使用自动导入功能后直接就去摸鱼，不然程序也会摸鱼的。  
  

2. 手动导入  
  当出现`手动导入数据集模式启动：`提示时 就说明配置文件中没有正确导入或者没有写入数据集地址，此时便要手动导入数据集。  
  手动导入数据集同ini配置导入，当数据集位于**本程序所在**的文件夹时使用相对路径，不在**本程序所在**的文件夹时使用绝对路径。 但是一次只能填写一个，例如：`C:\Users\xx\Desktop\xxx.tsv`或者`xxx.tsv`，输入完成后敲击回车即可。  
  当您输入完成最后一个数据集地址敲击回车后，可以输入`y`或者`Y`完成导入，至此您可以开始摸鱼了。  
  
#### **注意事项：**  
  1. GSE数据集请使用 ***.tsv*** 格式，并且保证其为 ***utf-8***编码格式  
  >编码格式查看可使用windows系统自带记事本软件打开数据集，右下角即为
  > 该tsv文件编码格式，若编码不是utf-8，可以将该文件另存为，在另存为时选择utf-8编码格式。
  
  2. GSE数据集解析分为两种模式  
  - 单集解析  
  单集解析即只导入一个数据集，当导入一个数据集时，**默认**读取**第一列**的数据， 所以请保证单集导入时，第一列数据为您的miRNA名称（hsa-miR-***格式）。  

  - 多集解析  
  多集解析即导入多个数据集，只要任意两个数据集之间的miRNA存在交集，便会进行查询，此时**默认**读取**第七列**的数据（hsa-miR-***格式）（大多数GSE数据集均第七列为miRNA名称）。  

  此外，请保证GSE数据集**第一行**为标题行，不然无法查询位于第一行的miRNA信息。  

- ### <span id='title2'>StarBase查询</span>  
  ***此步不需要操作，仅说明原理***  
>Starbase 是通过对高通量的 CLIP-Seq 实验数据和降解组实验数据来搜寻到 micorRNA 靶标，为探讨 microRNA 的靶标提供了各式各样的可视化界面，该数据库容括了丰富的 miRNA-ncRNA、miRNA-mRNA、RBP-RNA 和 RNA-RNA 的数据。

  当成功导入数据集后，程序便会第一个查询StarBase数据库，因为Ualcan数据库相关基因较少，且服务器位于境外，故本程序不采用Ualacn数据库为筛选数据库。  
  本程序在这一环节中会使用计算机爬虫技术，获取miRNA的差异表达和生存分析曲线的统计学p值，当两个p值**均小于等于**0.05，满足统计学意义时，便进入下一查询环节。
- ### <span id='title3'>Ualcan查询</span>  
  ***此步不需要操作，仅说明原理***  
>UALCAN是一个易于使用的，交互式的门户网站，可以执行对TCGA基因表达数据的深入分析。UALCAN使用TCGA level 3 RNA-seq和31种癌症的临床数据。

  当筛选完可用miRNA后，程序便会进一步查询Ualcan数据库，因为Ualcan数据库相关基因较少，可能存在只有前体没有成熟体或只有成熟体没有前体的情况，故本程序提供模糊搜索功能，将会根据miRNA名称，自动对其成熟体或前体再次查询，以达到全面获得数据的目的。  
  本程序在这一环节中会使用计算机爬虫技术，获取miRNA的差异表达和生存分析曲线的统计学p值，**无论**两个p值**是否均小于等于**0.05，都会进行记录并标明前体或成熟体。  
- ### <span id='title4'>中国知网&PubMed查询</span>  
  ***此步不需要操作，仅说明原理***  
  本程序会使用计算机爬虫技术对这两个数据库进行检索，获取该miRNA的相关文献，并对其进行词云解析。  

        'ACC': ['Adrenocortical Carcinoma', 'ACC', '肾上腺皮质癌'],  
        'BLCA': ['BLCA', 'Bladder Urothelial Carcinoma', '膀胱上皮癌', '膀胱癌', '膀胱尿路上皮癌'],
        'BRCA': ['BRCA', 'Breast Invasive Carcinoma', '乳腺癌', '乳腺浸润癌'],
        'CESC': ['CESC', 'Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma', '宫颈癌', '宫颈鳞状细胞癌和宫颈腺癌'],
        'CHOL': ['CHOL', 'Cholangiocarcinoma', '胆管癌', '胆管细胞癌'],
        'COAD': ['COAD', 'Colon Adenocarcinoma', '结肠腺癌', '大肠腺癌'],
        'DLBC': ['DLBC', 'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma', '淋巴肿瘤弥漫性大B细胞淋巴瘤'],
        'ESCA': ['ESCA', 'Esophageal Carcinoma', '食管癌'],
        'HNSC': ['HNSC', 'Head and Neck Squamous Cell Carcinoma', '头颈部鳞状细胞癌', '头颈部鳞癌'],
        'KICH': ['KICH', 'Kidney Chromophobe', '肾嫌色细胞瘤'],
        'KIRC': ['KIRC', 'Kidney Renal Clear Cell Carcinoma', '肾透明细胞癌'],
        'KIRP': ['KIRP', 'Kidney Renal Papillary Cell Carcinoma', '肾乳头细胞癌'],
        'LAML': ['LAML', 'Acute Myeloid Leukemia', '急性髓系白血病', '白血病'],
        'LGG': ['LGG', 'Brain Lower Grade Glioma', '脑低级别胶质瘤', '神经胶质瘤'],
        'LIHC': ['LIHC', 'Liver Hepatocellular Carcinoma', '肝细胞癌', '肝癌'],
        'LUAD': ['LUAD', 'Lung Adenocarcinoma', '肺腺癌'],
        'LUSC': ['LUSC', 'Lung Squamous Cell Carcinoma', '肺鳞状细胞癌'],
        'MESO': ['MESO', 'Mesothelioma', '间皮瘤'],
        'OV': ['OV', 'Ovarian Serous Cystadenocarcinoma', '卵巢浆液性囊腺癌', '卵巢癌'],
        'PAAD': ['PAAD', 'Pancreatic Adenocarcinoma', '胰腺癌'],
        'PCPG': ['PCPG', 'Pheochromocytoma and Paraganglioma', '嗜铬细胞瘤和副神经节瘤', '嗜铬细胞瘤', '副神经节瘤'],
        'PRAD': ['PRAD', 'Prostate Adenocarcinoma', '前列腺癌', '前列腺肿瘤'],
        'READ': ['READ', 'Rectum Adenocarcinoma', '直肠癌'],
        'SARC': ['SARC', 'Sarcoma', '肉瘤', '恶性毒瘤'],
        'SKCM': ['SKCM', 'Skin Cutaneous Melanoma', '皮肤黑色素瘤'],
        'STAD': ['STAD', 'Stomach Adenocarcinoma', '胃癌', '胃腺癌'],
        'TGCT': ['TGCT', 'Testicular Germ Cell Tumors', '睾丸生殖细胞瘤', '睾丸癌'],
        'THCA': ['THCA', 'Thyroid Carcinoma', '甲状腺癌', '甲状腺肿瘤'],
        'THYM': ['THYM', 'Thymoma', '胸腺癌', '胸腺肿瘤'],
        'UCEC': ['UCEC', 'Uterine Corpus Endometrial Carcinoma', '子宫体子宫内膜癌'],
        'UCS': ['UCS', 'Uterine Carcinosarcoma', '子宫癌肉瘤'],
        'UVM': ['UVM', 'Uveal Melanoma', '葡萄膜黑色素瘤']
***以上为关键词云***  
  我们在词云匹配使用了fuzzywuzzy模块，通过与关键词的匹配度计算degree。  
>  如果在进行关键词云匹配过程中程序终端有可能是因为您的计算机没有安装Visual C++ Build Tools for Visual Studio 2015 with Update 3。  
  <https://my.visualstudio.com/Downloads?q=Visual%20Studio%202015%20update%203&pgroup=>  
  请下载**DVD**版本，如果您的电脑储存不足，可以在下载后仅安装win10/win7 SDK.

-  当关键词云与相关论文摘要完全匹配时，degree为**Very High**，即有大概率这个miRNA有人做过相关癌种了，你做了大概率被导师骂一顿。  
-  当关键词云与相关论文摘要基本匹配时，degree为**High**，即有一定概率这个miRNA有人做过相关癌种了。  
-  当关键词云与相关论文摘要基本不匹配时（可能有其他癌种有相关研究），degree为**Low**，即有小概率这个miRNA有人做过相关癌种了。  
-  当没有相关论文时，degree为**None**，即基本没人做过这个基因的任何研究。  

  程序会自动对miRNA进行分类，根据米娜桑在配置文件中设置的`model`值进行下一步搜索。
- ### <span id='title6'>mirWalk查询</span>  
  ***此步不需要操作，仅说明原理***  
>  miRWalk（ <http://mirwalk.umm.uni-heidelberg.de/> ）是一个综合性的miRNA靶基因数据库，收录了Human、Mouse、Rat、Dog、cow等多个物种的miRNA靶基因信息，不仅仅记录了基因全长序列上的miRNA结合位点，也会将其与已有的12个miRNA靶标预测程序（DIANA-microTv4.0 , DIANA-microT-CDS , miRanda-rel2010 , mirBridge , miRDB4.0 , miRmap , miRNAMap, doRiNA i.e.,PicTar2, PITA RNA22v2 , RNAhybrid2.1 and Targetscan6.2 ）的预测结合信息集合进行结合关联。数据库一直在更新收录新的资料，第一版是在2011年发布，之后在2015年发布V2 版本，并登上了Nature methods杂志，目前更新到V3 版本。

  根据米娜桑填写的`model`值，程序会查询mirWalk数据库，程序通过pandas模块对下载到的csv文件进行解析，获取靶蛋白`bindingp` = 1的靶蛋白。  
- ### <span id='title5'>miRDB查询</span>  
  ***此步不需要操作，仅说明原理***  
>  miRDB（http://mirdb.org）是用于miRNA目标预测和功能注释的在线数据库。miRDB中的所有目标都由生物信息学工具MirTarget预测，该工具是通过分析来自高通量测序实验的数千个miRNA-目标相互作用而开发的。与miRNA结合和靶点下调相关的共同特征已被确定并用于通过机器学习方法预测miRNA靶点。

  根据米娜桑填写的`model`值，程序会查询miRDB数据库，程序通过xpath模块对页面进行解析，获取靶蛋白`Target Score` > 80的靶蛋白。
- ### <span id='title6'>TargetScan查询</span>  
  ***此步不需要操作，仅说明原理***  
>  TargetScan （http://www.targetscan.org/vert_71/） 是一个miRNA 靶基因预测的网站, 包括了 人， 小鼠，果蝇 ， 线虫， 斑马鱼 共5个物种的miRNA 靶基因结果。

  根据米娜桑填写的`model`值，程序会查询TargetScan数据库，程序通过panads模块对下载到的xlsx文件进行解析，获取靶蛋白`Total context++ score` < -0.5 的靶蛋白。  
- ### <span id='title7'>mirDIP查询</span>  
  ***此步不需要操作，仅说明原理***  
>  mirDIP（<http://ophid.utoronto.ca/mirDIP/index.jsp>）集成了30个来源数据库中human相关的靶基因信息，是最全面的人类miRNA靶基因数据库。

  根据米娜桑填写的`model`值，程序会查询mirDIP数据库，程序通过对下载到的tsv文件进行解析，获取靶蛋白`very high`的靶蛋白。
- ### <span id='title8'>Tarbase查询</span>  
  ***此步不需要操作，仅说明原理***  
>  TarBase（<https://dianalab.e-ce.uth.gr/html/diana/web/index.php?r=tarbasev8>）收录各种实验验证过的miR-Target数据（只要是人和小鼠的靶基因信息），TarBase将实验证据分为low和high两类，low代表的是传统的实验手段，可靠性相对于高通量测序的分析结果更高一点，我们可以筛选low方法支持的miRNA靶基因信息，得到高质量的miRNA候选靶基因集。TarBase只提供在线检索，可以输入miRNA名称和/或基因名称，miRNA名称为miRBase数据库格式,基因名称支持gene symbol(基因名)和ensembl gene ID。

  根据米娜桑填写的`model`值，程序会查询Tarbase数据库，程序通过xpath模块对页面进行解析，获取靶蛋白`Score` > 5 的靶蛋白。
