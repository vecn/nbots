LIB_NAME = jni_vcn_bots
VCN_LIB = -L ../../shared_libs/
VCN_INCLUDE = -I ../../include/
JDK_INCLUDE = -I /usr/local/java/jdk1.8.0_45/include/ -I /usr/local/java/jdk1.8.0_45/include/linux/

all:
	gcc -std=c99 -c -fPIC ${JDK_INCLUDE} ${VCN_INCLUDE} jni_geometric_bot.c -o ${LIB_NAME}.o
	gcc -shared ${VCN_LIB} -o lib${LIB_NAME}.so ${LIB_NAME}.o -lvcn_bots
	rm -f ${LIB_NAME}.o