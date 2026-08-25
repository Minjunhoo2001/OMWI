import pandas as pd
import os
import sys

def extract_metaphlan_levels(input_file, output_dir):
    """
    从MetaPhlAn合并结果中提取不同层级的分类丰度
    参数:
        input_file: 输入的TSV文件路径
        output_dir: 输出目录路径
    """
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取文件，跳过注释行
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    # MetaPhlAn merged tables use `clade_name`; some older exports use `ID`.
    # Fall back to the first column so the extractor remains compatible with
    # both formats while preserving the original taxon column in its output.
    taxon_col = next(
        (name for name in ('clade_name', '#clade_name', 'ID') if name in df.columns),
        df.columns[0] if len(df.columns) else None,
    )
    if taxon_col is None:
        raise ValueError('输入表为空或缺少分类名称列')

    # 计算分类层级深度
    df['level'] = df[taxon_col].astype(str).str.count(r'\|') + 1
    
    # 定义层级标签
    level_labels = {
        1: 'kingdom',
        2: 'phylum',
        3: 'class',
        4: 'order',
        5: 'family',
        6: 'genus',
        7: 'species',
        8: 'SGB'
    }
    
    # 按层级分组并保存
    for level_num, level_name in level_labels.items():
        # 提取当前层级数据
        level_df = df[df['level'] == level_num].copy()
        
        if not level_df.empty:
            # 移除临时添加的层级列
            level_df.drop(columns=['level'], inplace=True)
            
            # 生成输出文件名
            base_name = os.path.basename(input_file).split('.')[0]
            output_file = os.path.join(output_dir, f"{base_name}_{level_name}.tsv")
            
            # 保存结果
            level_df.to_csv(output_file, sep='\t', index=False)
            print(f"已创建: {output_file} ({len(level_df)}行)")
        else:
            print(f"警告: {level_name}层级未找到数据")

def main():
    if len(sys.argv) != 3:
        print("用法: python metaphlan_level_extractor.py <输入文件.tsv> <输出目录>")
        print("示例: python metaphlan_level_extractor.py taxonomy.tsv level_results")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.isfile(input_file):
        print(f"错误: 输入文件不存在 {input_file}")
        sys.exit(1)
    
    extract_metaphlan_levels(input_file, output_dir)
    print("处理完成！")

if __name__ == "__main__":
    main()
