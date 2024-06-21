#! /bin/bash

# 定义obabel路径
obabel=/pub/anaconda3/bin/obabel

# 加载环境变量
source /etc/profile

# 转换SDF文件为PDB文件
$obabel -isdf pNPR.sdf -opdb -O ligand.pdb

# 准备受体和配体文件
prepare_receptor -r rhapaa.pdb -o protein.pdbqt -A hydrogens
prepare_ligand -l ligand.pdb -o ligand.pdbqt

# 初始化循环变量
i=1
index=1
end=20

# 创建存放PDB文件的目录
mkdir pdb

# 主循环，运行50次对接
while (( i <= 50 )); do
    # 执行对接并保存输出
    vina --config config.txt --out dock_$i.pdbqt

    # 创建子目录并移动对接结果
    mkdir dock_$i
    cp dock_$i.pdbqt dock_$i/

    # 进入子目录并拆分对接结果
    cd dock_$i
    vina_split --input dock_$i.pdbqt
    cd ..

    # 内部循环，处理每个对接结果
    while (( index <= end )); do
        t=$(printf "%02d" $((index + 20 - end)))
        echo $obabel -ipdbqt dock_$i/dock_${i}_ligand_$t.pdbqt -opdb -O pdb/${index}.pdb
        $obabel -i pdbqt dock_$i/dock_${i}_ligand_$t.pdbqt -o pdb -O pdb/${index}.pdb
        index=$((index + 1))
    done

    # 更新结束索引和主循环计数器
    end=$((end + 20))
    i=$((i + 1))
done
