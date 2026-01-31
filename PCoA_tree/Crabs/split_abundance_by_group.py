import pandas as pd

# 1. 读取文件
# 丰度文件：第一行是样本ID，第一列是物种分类
abundance_df = pd.read_csv("family_out1.txt", sep="\t", index_col=0)  # 以第一列为索引（物种）
# 分组文件：读取后构建“样本ID→分组”的字典
metadata_df = pd.read_csv("metadata1.txt", sep="\t")
sample_group_dict = dict(zip(metadata_df["sample-id"], metadata_df["gutregion"]))

# 2. 筛选每个分组对应的样本列
# 获取所有分组（H/M/F）
groups = ["H", "M", "F"]
for group in groups:
    # 找到当前分组对应的所有样本ID
    group_samples = [sample for sample, g in sample_group_dict.items() if g == group]
    # 筛选丰度文件中对应的列（保留物种列+分组样本列）
    group_abundance = abundance_df[group_samples]
    # 3. 保存为新文件（保留原始格式）
    output_file = f"family_out1_{group}.txt"
    group_abundance.to_csv(output_file, sep="\t")
    print(f"已生成分组文件：{output_file}")

print("所有分组文件已生成完成！")
