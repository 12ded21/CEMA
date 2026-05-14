#include "utils.h"

extern int K;
extern int *degree;
extern int *dedges;

std::string integeter_to_str(ll number){
    std::vector<int> sequence;
	if(number == 0) 
        sequence.pb(0);
	while(number > 0) {
		sequence.pb(number%1000);
		number /= 1000;
	}

	char buf[5];
	std::string res;
	for(int i = sequence.size();i > 0;i --) {
		if(i == sequence.size()) 
            sprintf(buf, "%u", sequence[i-1]);
		else 
            sprintf(buf, ",%03u", sequence[i-1]);
		res += std::string(buf);
	}
	return res;
}

FILE* open_file(const char* file_name, const char* mode) {
	FILE* f = fopen(file_name, mode);
	if (f == nullptr) {
		printf("Can not open file: %s\n", file_name);
		exit(1);
	}
	return f;
}

void fread_wall(void* __ptr, size_t __size, size_t __n, FILE* __stream) {
	int sz = fread(__ptr, __size, __n, __stream);
	if (sz != __n) {
		fprintf(stderr, "Error in reading\n");
	}
}

int func(int l, int q){
    if(l > q)
        return 0;
    else{
        int a = std::ceil(((double)q) / ((double)(l - 1)));
        int b = (l - 1) * a - q;
        int c = (q - b * (a - 1)) / a;
        int d = (a - 1)*(a - 2) / 2;
        int e = a * (a - 1) / 2;
        return b * d + c * e;
    }
}

int f1(int budget, int clq_size){
	return std::min((int)ceil(log2(budget + 1)), clq_size * (clq_size - 1) / 2);
}

int f2(int budget, int clq_size, int del, bool flag){
	if(del < 0){
		del = std::min((int)ceil(log2(budget + 1)), clq_size * (clq_size - 1) / 2);
	}
	else{
		if(flag){
			del = std::min(clq_size * ( clq_size - 1 ) / 2, del + 1);
		}
		else{
			del = std::max(1, del - 1);
		}
	}
	return del;
}

int f3(int budget, int clq_size, int del, bool flag){
	if(del < 0){
		del = 1;
	}
	else{
		if(flag){
			del = 1;
		}
		else{
			del = std::min((int)(del + 1), clq_size * (clq_size - 1) / 2);
		}
	}
	return std::min(del, clq_size * (clq_size - 1) / 2);
}
