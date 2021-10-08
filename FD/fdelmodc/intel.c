extern  int __intel_cpu_indicator;
extern  void __intel_cpu_indicator_init(void)  ;
extern  void __intel_new_proc_init_P(void)  ;
extern  void __intel_new_proc_init(void)  ;
__intel_cpu_indicator_init();
__intel_cpu_indicator = 2048; /* 2048 for EMT64, 1024 for pentium-m, 512 for xeon...*/
//__intel_cpu_indicator = -512; /* 2048 for EMT64, 1024 for pentium-m, 512 for xeon...*/

void __intel_cpu_indicator_init(void){
    return;
}

void __intel_new_proc_init_P(void){
    return;
}
void __intel_new_proc_init(void){
    return;
}

