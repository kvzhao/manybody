RESFILE=$1
python3 plot_observable_vs_size.py --observ=E --data_dir=${RESFILE} --out_name=figures/Energy
python3 plot_observable_vs_size.py --observ=absM --data_dir=${RESFILE} --out_name=figures/AbsoluteMagnetization
python3 plot_observable_vs_size.py --observ=M2 --data_dir=${RESFILE} --out_name=figures/Magnetization
python3 plot_observable_vs_size.py --observ=Cv --data_dir=${RESFILE} --out_name=figures/HeatCapacity
python3 plot_observable_vs_size.py --observ=U2 --data_dir=${RESFILE} --out_name=figures/BinderRatio