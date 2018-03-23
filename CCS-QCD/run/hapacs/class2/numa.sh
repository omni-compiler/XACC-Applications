LOCAL_RANK=$MV2_COMM_WORLD_LOCAL_RANK
SOCKET=$(expr $LOCAL_RANK / 2)
numactl --cpunodebind=$SOCKET --localalloc $@
