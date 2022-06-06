#' getlcai
#'
#' Using the mRNA expression profile data and sample grouping information data,
#' you can calculate the Lung Cancer Aggressive Index (LCAI).
#'
#' @param exp a expression matrix (\code{data.frame} or \code{file}).
#' The columns of the matrix are samples and the rows are the unique genes.
#'
#' @param pheno a sample grouping information matrix (\code{data.frame} or \code{file}).
#' There are two columns in the matrix. The first column is the sample ID and
#' the second column is the grouping information.
#'
#' @param control a character vector. The group you want to set as control in
#' the group information.
#'
#' @param experimental a character vector. The group you want to set as
#' experimental in the group information.
#'
#' @param type the type of experiment to obtain mRNA expression profile data.
#' Use \code{"Array"} or \code{"RNA-seq"}.
#'
#' @param plotPCA a boolean defaulting to FALSE as to whether or not to plot 3D
#' PCA graph.Use \code{TRUE} or \code{FALSE}.
#'
#' @description Calculate the Lung Cancer Aggressive Index (LCAI):
#' This requires the following inputs:
#'
#' The mRNA_expression_matrix (exp):
#'
#'                   GSM5342379 GSM5342380 GSM5342381 GSM5342382 GSM5342383 ...
#'
#'      LOC100130938  0.9945589  0.9313970  2.1251374  1.7672035  0.9419411
#'
#'      LOC100507487  1.1168837  1.6486066  1.2722476  1.0819411  1.4785424
#'
#'      LOC145474     0.8186053  0.7759738  0.8175747  0.7968981  0.9242794
#'
#'      LOC102723721  1.7684608  0.7185605  0.6933824  0.6953136  2.2042653
#'
#'      KIAA0040      0.9653091  1.2400747  1.5995581  1.8094245  2.9599967
#'
#'      LOC101927686  1.2996664  0.9895734  1.8646964  1.5624728  2.1886328
#'
#'      ...
#'
#'
#' The sample grouping information matrix (pheno), contains sample names and group names:
#'
#'           sample        group
#'
#'      1 GSM5342379      Control
#'
#'      2 GSM5342380      Control
#'
#'      3 GSM5342381      Control
#'
#'      4 GSM5342382 experimental
#'
#'      5 GSM5342383 experimental
#'
#'      6 GSM5342384 experimental
#'
#'      ...
#'
#' @return a list: the lcai value; mergerd expression matrix; group information;the first three principal components;
#' @export
#'
#' @examples 
#' data(exp_example)
#' data(pheno_example)
#' outlist = getlcai(exp_example,pheno_example,"Control","experimental",type='Array',plotPCA = TRUE)
#'
#'
#' @import dplyr
#' @import rgl
#' @importFrom methods is
#' @importFrom stats prcomp aggregate
#' @importFrom utils read.table

getlcai <- function(exp,
                    pheno,
                    control,
                    experimental,
                    type='Array',
                    plotPCA = FALSE){


  # exp = "../../GSE175601/GSE175601_unique_exp.txt"
  # pheno = "../../GSE175601/GSE175601_phe.txt"
  # control = "Control"
  # experimental = "experimental"
  # plotPCA = TRUE
  # 声名函数内变量
  group = NULL
  color = NULL
  pheno_pc = NULL
  fig = NULL

  message("[===========================]")
  message("[###### getlcai START ######]")
  message("-----------------------------")

  

  #### 0 从用户获取输入数据 ====
  message("\n[INFO] Reading the input data. Please ensure that the input data format is consistent with the requirements.")

  # 读取输入表达矩阵 data.frame或文本文件 第一列为基因名,其余列为样本名
  if (Reduce("|",is(exp) == "character")) {
    # 判断input文件路径是否正确
    if(file.exists(exp) == FALSE){
      stop("Reading error: Your 'exp' for loading does not exists, please assign a correct file path.")
    }else{
      exp_input = read.table(exp,
                             stringsAsFactors = F,
                             sep = "\t",
                             header = T,
                             row.names = 1)
    }
  }else if (Reduce("|",is(exp) %in% c("data.frame"))){
    exp_input = exp
  }else if (Reduce("|",is(exp) %in% c("matrix"))){
    exp_input = as.data.frame(exp)
  }else{
    stop("Reading error: 'exp' isn't recognized as a data.frame or filename.")
  }

  # 读取输入分组信息 data.frame或文本文件
  if (Reduce("|",is(pheno) == "character")) {
    # 判断input文件路径是否正确
    if(file.exists(pheno) == FALSE){
      stop("Reading error: Your 'pheno' for loading does not exists, please assign a correct file path.")
    }else{
      pheno_input = read.table(pheno,
                               stringsAsFactors = F,
                               sep = "\t",
                               header = T)
    }
  }else if (Reduce("|",is(pheno) %in% c("data.frame"))){
    pheno_input = pheno
  }else if (Reduce("|",is(pheno) %in% c("matrix"))){
    pheno_input = as.data.frame(pheno)
  }else{
    stop("Reading error: 'pheno' isn't recognized as a data.frame or filename.")
  }

  message("[INFO] Reading Success")

  #### 1 处理输入矩阵,只保留mapping到的基因的表达量 ====
  # 判断输入矩阵测序方法：Array、RNA-seq
  # 重点：这里涉及到一个不同技术平台数据合并的问题，目前只是单纯的合并，会有问题，后续需要更好的合并标准化方法
  exp_input = as.data.frame(lapply(exp_input,as.numeric),row.names = rownames(exp_input))

  if (type == "RNA-seq") {
    exp_input = log(exp_input+1)
    message("[INFO] Data type is RNA-seq.The data is transformed by natural logarithm.")
  }else if (type == "Array") {
    message("[INFO] Data type is Array.")
  }else{
    stop("Error: Data type setting error! Please reset it by 'RNA SEQ' or 'Array'.")
  }

  # 输入矩阵按照背景矩阵基因id进行mapping
  # exp_input = exp_input[as.factor(rownames(exp_input)) %in% id2entrez$symbol, ]             # 判断并保留exp在背景数据的基因所在行
  # exp_input$ids = id2entrez[match(as.factor(rownames(exp_input)), id2entrez$symbol), 2]     # 映射 基因id 到 Entrez id

  # 将基因名变为小写，再进行匹配
  exp_input = exp_input[tolower(as.factor(rownames(exp_input))) %in% tolower(id2entrez$symbol), ]             # 判断并保留exp在背景数据的基因所在行
  exp_input$ids = id2entrez[match(tolower(as.factor(rownames(exp_input))), tolower(id2entrez$symbol)), 2]     # 映射 基因id 到 Entrez id

  # symbol转化为Entrez id后可能会出现重复，加一个去重
  exp_input = aggregate(.~ids,exp_input,mean)   # 去重，取平均值
  rownames(exp_input) = exp_input$ids
  exp_input = exp_input[,-1]                    # 行名改为 Entrez id



  #### 2 输入矩阵与背景矩阵合并 ====
  exp_input$id = rownames(exp_input)           # 添加一列id,后续基于id列合并
  exp_back$id = rownames(exp_back)
  exp_join = dplyr::inner_join(exp_back, exp_input, by = 'id')      # 按照id列,将exp_input中对应的数据添加到dat末尾列后,inner_join取交集

  # 输出交集基因数
  gene_left_count = nrow(exp_back)
  gene_right_count = nrow(exp_join)
  gene_percent = paste(sprintf("%.2f", (gene_right_count / gene_left_count) * 100), "%", sep = "")              # 计算百分比
  gene_info = data.frame(GEO_genes = gene_left_count,
                         merged_genes = gene_right_count,
                         proportion = gene_percent)

  # 输出信息 背景数据基因数, 合并后基因数, 占比
  message(
    "\n[INFO] The number of GEO background genes is ",
    gene_left_count,
    ", the number of genes after merged is ",
    gene_right_count,
    ", account for ",
    gene_percent,
    "\n"
  )
  print(gene_info,quote=F)

  # 完善exp_join格式
  rownames(exp_join) = exp_join$id              # 修改行名为id
  exp_join = dplyr::select(exp_join, -id)       # 删去id那一列

  # if (type == "RNA-seq") {
  #   exp_join = log(exp_join+1)
  #   message("[INFO] Data type is RNA-seq.The data is transformed by natural logarithm.")
  # }

  
  
  #### 3 样本分组信息合并 ====
  colnames(pheno_input) = c("sample", "group")

  # 判断输入文件中分组名是否和背景数据分组名相同
  if (any(pheno_input$group %in% c("NSCLC", "SCLC", "NLE"))) {
    stop("Error: Please do not set 'NSCLC,SCLC or NLE' as group names.")
  } else{
    pheno = rbind(pheno_back, pheno_input)
  }

  
  
  #### 4 PCA 主成分分析 ====

  # 标准化,这里使用R自带的scale函数 其他方法待定
  # exp_join = data.frame(scale(exp_join), check.names = F)
  pc = prcomp(t(exp_join),cente=TRUE,scale.=TRUE)
  pcdata = data.frame(pc$x[, 1:3])

  # 循环计算每个分组中心坐标(平均值)
  PC_group_mean = c("group", "PC1", "PC2", "PC3")
  phe_unique = dplyr::count(x = pheno, group)

  for (g in phe_unique$group) {
    dat = pcdata[pheno[pheno$group == g, ]$sample, ]
    PC_group_mean = rbind(PC_group_mean, c(g, mean(dat$PC1), mean(dat$PC2), mean(dat$PC3)))
  }

  PC_group_mean = data.frame(PC_group_mean,
                             row.names = NULL,
                             stringsAsFactors = F)
  colnames(PC_group_mean) = PC_group_mean[1, ]

  # 每个分组的中心点
  group_center = PC_group_mean[-1, ]


  #### 5 肺癌计算进展量化指数 ====
  NSCLC_center = group_center[group_center$group == "NSCLC", ]
  SCLC_center = group_center[group_center$group == "SCLC", ]
  group_former_center = group_center[group_center$group == control, ]
  group_later_center = group_center[group_center$group == experimental, ]

  ## 进展向量: NSCLC_center ===> SCLC_center
  vector_back = get_vector(NSCLC_center,SCLC_center)

  ## 处理前向量(对照): NSCLC_center ===> group_former_center
  vector_v = get_vector(NSCLC_center,group_former_center)
  # 向量欧氏距离:
  Lv = sqrt(vector_v[1]^ 2+vector_v[2]^ 2+vector_v[3]^ 2)
  # 投影v
  cosv = get_cos(vector_v, vector_back)
  Lv1 = Lv * cosv

  ## 处理后向量(实验): NSCLC_center ===> group_later_center
  vector_d = get_vector(NSCLC_center,group_later_center)
  # 向量欧氏距离:
  Ld = sqrt(vector_d[1]^ 2+vector_d[2]^ 2+vector_d[3]^ 2)
  # 投影d
  cosd = get_cos(vector_d, vector_back)
  Ld1 = Ld * cosd

  ## LCAI
  lcai = Ld1 - Lv1


  
  #### 6 绘制3D PCA图 ====

  if(plotPCA){
    # 颜色矩阵，后续再完善下
    colorname = c("blue","hotpink","green","yellow","red")
    group_count = length(table(pheno$group))
    pheno$color = factor(pheno$group, labels = colorname[1:group_count])

    # 画布大小
    rgl::open3d(windowRect=c(100,100,700,700))
    rgl::plot3d(pcdata,
                col = pheno$color,
                type = 'p',
                size = 8)
    # 绘制直线叠加到上一个图
    rgl::plot3d(x = c(group_former_center$PC1,group_later_center$PC1),
                y = c(group_former_center$PC2,group_later_center$PC2),
                z = c(group_former_center$PC3,group_later_center$PC3),
                col = c("black","black"),add = T,type = 'p',size = 8)

    rgl::plot3d(x = c(group_former_center$PC1,group_later_center$PC1),
                y = c(group_former_center$PC2,group_later_center$PC2),
                z = c(group_former_center$PC3,group_later_center$PC3),
                col = "black",add = T,type = 'l',lwd = 2)

    # pheno_pc = cbind(pheno,pcdata)
    # fig = plotly::plot_ly(pheno_pc,x=~PC1,y=~PC2,z=~PC3,color = ~group,colors = ~c("blue","hotpink","green","yellow","red"),size = 2)
    # fig = plotly::add_markers(fig)

    # 输出信息
    phe_unique = dplyr::count(x = pheno, group, color)
    message("\n[INFO] PCA graph plotting complete! Colors and groups information are as follows\n")
    print(phe_unique)
  }


  
  #### 7 输出LCAI值及其他信息 ====

  # 判断正负并输出lcai到屏幕
  if (lcai > 0) {
    lcai = paste("+", lcai, sep = "")
    message(
      "\n[INFO] The lung cancer aggressive index (LCAI) is ",
      lcai,
      ", group: ",
      control,
      " ---> ",
      experimental
      # ", showing the cancer promoting effect."
    )
  } else {
    message(
      "\n[INFO] The lung cancer aggressive index (LCAI) is ",
      lcai,
      " [",
      control,
      " ---> ",
      experimental,
      "]"

      # ", showing the cancer inhibition effect."
    )
  }

  message("\n")
  message("[===========================]")
  message("[######  getlcai END  ######]")
  message("-----------------------------")

  # 输出 lcai, 分组数据,前3个PC数据
  out = list(lcai=lcai,exp=exp_join,pheno=pheno,pcdata=pcdata)
  # out = list(lcai=lcai,exp=exp_join,pheno_pc=pheno_pc,fig=fig)
  return(out)
}







#' 计算向量坐标
#'
#' @param center1 a numeric vector
#' @param center2 a numeric vector
#' @noRd
get_vector = function(center1, center2){
  center1 = as.numeric(center1[, -1])
  center2 = as.numeric(center2[, -1])
  vector = c(center2[1] - center1[1],center2[2] - center1[2],center2[3] - center1[3])
  return(vector)
}


#' 计算向量夹角余弦
#'
#' @param vector1 a numeric vector
#' @param vector2 a numeric vector
#' @noRd
get_cos = function(vector1, vector2) {
  dot1 = vector1[1] * vector2[1] + vector1[2] * vector2[2] + vector1[3] * vector2[3]
  dist1 = sqrt(vector1[1] ^ 2 + vector1[2] ^ 2 + vector1[3] ^ 2)
  dist2 = sqrt(vector2[1] ^ 2 + vector2[2] ^ 2 + vector2[3] ^ 2)
  cos = dot1 / (dist1 * dist2)
  return(cos)
}
