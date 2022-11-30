#splitfa $work_dir/PRC34mod.psmcfa > $work_dir/split.psmcfa
#psmc -N30 -t15 -r5 -p "64*1" -o $work_dir/PRC34.psmc $work_dir/PRC34mod.psmcfa
#seq 100 | xargs -n 1 -p 32 -i echo psmc -N30 -t15 -r5 -b -p "64*1" -o $work_dir/round-{}.psmc $work_dir/split.psmcfa | sh
#cat $work_dir/PRC34.psmc $work_dir/round-*.psmc > $work_dir/combined.psmc
#rm -vrf $work_dir/split.psmcfa $work_dir/PRC34.psmc
#mv $work_dir/combined.psmc $work_dir/PRC34.psmc
#mv $work_dir/PRC34.psmc $work_dir/psmc
#mv $work_dir/PRC34mod.psmcfa $work_dir/psmcfa
psmc -N30 -t15 -r5 -b -p "64*1" -o $2 $1
