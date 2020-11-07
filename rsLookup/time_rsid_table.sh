#!/bin/bash
i=0
real_time=()
user_time=()
sys_time=()

array_sum() {
	sum_total=0
	arr=("$@")
	for i in ${arr[@]}; do
		sum_total=$(bc <<< "$i + $sum_total")
	done
	echo $sum_total
}

convert_time() {
	time_s=${1%?}
	IFS='m' read -r -a t_a <<< "$time_s"
	min=$(bc <<< 60*${t_a[0]})
	tt=$(bc <<< "$min + ${t_a[1]}")
	echo $tt
}

if [ $1 == "rsid" ]
then
while IFS= read -r rsid
do
# { tte=$( { time ./genome_hash.exe get_rsid rsid.dht $rsid 1>&3- 2>&4-; } 2>&1 ); } 3>&1 4>&2
./genome_hash.exe get_rsid rsid.dht $rsid
: <<'END'
t_arr=()
while read -r t_str
do
read -r -a t_num <<< "$t_str"
t_arr+=("${t_num[1]}")
done <<< "$tte"
real_t=${t_arr[1]}
user_t=${t_arr[2]}
sys_t=${t_arr[3]}
real_time+=($(convert_time ${real_t}))
user_time+=($(convert_time ${user_t}))
sys_time+=($(convert_time ${sys_t}))
i=$((i+1))
END
done < rsid_lookup_tests.in
elif [ $1 == "snp_data" ]
then
while IFS= read -r line
do
chromosome=$(echo "$line" | cut -f1 -d$'\t')
position=$(echo "$line" | cut -f2 -d$'\t')
allele=$(echo "$line" | cut -f3 -d$'\t')
./genome_hash.exe get_cpa snp.dht $chromosome $position $allele
done < snp_tests.in
elif [ $1 == "snp_pointer" ]
then
while IFS= read -r line
do
chromosome=$(echo "$line" | cut -f1 -d$'\t')
position=$(echo "$line" | cut -f2 -d$'\t')
allele=$(echo "$line" | cut -f3 -d$'\t')
./genome_hash.exe get_cpa_pointer ../../../snp150Common.txt snp_pointer.dht $chromosome $position $allele
done < snp_tests.in
fi

exit

real_sum=$(array_sum "${real_time[@]}")
user_sum=$(array_sum "${user_time[@]}")
sys_sum=$(array_sum "${sys_time[@]}")

real_avg=$(bc -l <<< ${real_sum}/${i})
user_avg=$(bc -l <<< ${user_sum}/${i})
sys_avg=$(bc -l <<< ${sys_sum}/${i})

echo "Out of $i tests run:"
echo "Average real time: $real_avg"
echo "Average user time: $user_avg"
echo "Average sys time: $sys_avg"
