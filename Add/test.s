	.section	__TEXT,__text,regular,pure_instructions
	.globl	_TAB
	.align	4, 0x90
_TAB:
Leh_func_begin1:
	pushq	%rbp
Ltmp0:
	movq	%rsp, %rbp
Ltmp1:
	leaq	-56(%rbp), %rax
	movabsq	$44, %rcx
	addq	%rcx, %rax
	movq	%rax, -16(%rbp)
	movq	-16(%rbp), %rax
	movq	%rax, -8(%rbp)
	movq	-8(%rbp), %rax
	popq	%rbp
	ret
Leh_func_end1:

	.globl	_Victim
	.align	4, 0x90
_Victim:
Leh_func_begin2:
	pushq	%rbp
Ltmp2:
	movq	%rsp, %rbp
Ltmp3:
	subq	$16, %rsp
Ltmp4:
	callq	_TAB
	movq	%rax, %rcx
	movq	%rcx, -8(%rbp)
	movq	-8(%rbp), %rcx
	movl	$42, (%rcx)
	movq	-8(%rbp), %rcx
	movl	(%rcx), %ecx
	movq	___stderrp@GOTPCREL(%rip), %rsi
	movq	(%rsi), %rsi
	xorb	%dil, %dil
	leaq	L_.str(%rip), %r8
	movb	%dil, -9(%rbp)
	movq	%rsi, %rdi
	movq	%r8, %rsi
	movl	%ecx, %edx
	movb	-9(%rbp), %cl
	movb	%cl, %al
	callq	_fprintf
	addq	$16, %rsp
	popq	%rbp
	ret
Leh_func_end2:

	.section	__TEXT,__literal8,8byte_literals
	.align	3
LCPI3_0:
	.quad	4609434218613702656
LCPI3_1:
	.quad	4985484787500187648
LCPI3_2:
	.quad	4890909195324358656
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_main
	.align	4, 0x90
_main:
Leh_func_begin3:
	pushq	%rbp
Ltmp5:
	movq	%rsp, %rbp
Ltmp6:
	subq	$224, %rsp
Ltmp7:
	movl	%edi, %eax
	movl	%eax, -4(%rbp)
	movq	%rsi, -16(%rbp)
	movq	$2500, -72(%rbp)
	movq	$100000000, -64(%rbp)
	movq	$32, -56(%rbp)
	movq	-56(%rbp), %rax
	movabsq	$4, %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	movq	%rcx, -152(%rbp)
	callq	_calloc
	movq	%rax, -96(%rbp)
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	movq	-152(%rbp), %rsi
	callq	_calloc
	movq	%rax, -104(%rbp)
	movq	-56(%rbp), %rax
	movq	%rax, %rdi
	movq	-152(%rbp), %rsi
	callq	_calloc
	movq	%rax, -112(%rbp)
	movq	$0, -40(%rbp)
	jmp	LBB3_2
LBB3_1:
	callq	_rand
	movl	%eax, %ecx
	cvtsi2ss	%ecx, %xmm0
	movabsq	$2147483648, %rcx
	cvtsi2ssq	%rcx, %xmm1
	divss	%xmm1, %xmm0
	movq	-96(%rbp), %rcx
	movq	-40(%rbp), %rdx
	movss	%xmm0, (%rcx,%rdx,4)
	movss	%xmm1, -156(%rbp)
	callq	_rand
	movl	%eax, %ecx
	cvtsi2ss	%ecx, %xmm0
	movss	-156(%rbp), %xmm1
	divss	%xmm1, %xmm0
	movq	-104(%rbp), %rcx
	movq	-40(%rbp), %rdx
	movss	%xmm0, (%rcx,%rdx,4)
	callq	_rand
	movl	%eax, %ecx
	cvtsi2ss	%ecx, %xmm0
	movss	-156(%rbp), %xmm1
	divss	%xmm1, %xmm0
	movq	-112(%rbp), %rcx
	movq	-40(%rbp), %rdx
	movss	%xmm0, (%rcx,%rdx,4)
	movq	-40(%rbp), %rcx
	movabsq	$1, %rdx
	addq	%rdx, %rcx
	movq	%rcx, -40(%rbp)
LBB3_2:
	movq	-40(%rbp), %rax
	movq	-56(%rbp), %rcx
	cmpq	%rcx, %rax
	jb	LBB3_1
	callq	_wallclock_time
	movsd	%xmm0, -120(%rbp)
	movq	$0, -32(%rbp)
	jmp	LBB3_8
LBB3_4:
	movq	$0, -48(%rbp)
	jmp	LBB3_6
LBB3_5:
	movq	-104(%rbp), %rax
	movq	-48(%rbp), %rcx
	movss	(%rax,%rcx,4), %xmm0
	cvtss2sd	%xmm0, %xmm0
	movsd	LCPI3_0(%rip), %xmm1
	callq	_pow
	cvtsd2ss	%xmm0, %xmm0
	movq	-112(%rbp), %rax
	movq	-48(%rbp), %rcx
	movss	%xmm0, (%rax,%rcx,4)
	movq	-48(%rbp), %rax
	movabsq	$1, %rcx
	addq	%rcx, %rax
	movq	%rax, -48(%rbp)
LBB3_6:
	movq	-48(%rbp), %rax
	movq	-56(%rbp), %rcx
	cmpq	%rcx, %rax
	jb	LBB3_5
	movq	-32(%rbp), %rax
	movabsq	$1, %rcx
	addq	%rcx, %rax
	movq	%rax, -32(%rbp)
LBB3_8:
	movq	-32(%rbp), %rax
	movq	-64(%rbp), %rcx
	cmpq	%rcx, %rax
	jb	LBB3_4
	callq	_wallclock_time
	movsd	%xmm0, -128(%rbp)
	movq	-64(%rbp), %rax
	movq	-56(%rbp), %rcx
	movq	%rcx, %rdx
	shrq	$32, %rdx
	movabsq	$4985484787499139072, %rsi
	leaq	(%rdx,%rsi), %rdx
	movd	%rdx, %xmm0
	movsd	LCPI3_1(%rip), %xmm1
	subsd	%xmm1, %xmm0
	movl	%ecx, %ecx
	movabsq	$4841369599423283200, %rdx
	leaq	(%rcx,%rdx), %rcx
	movd	%rcx, %xmm2
	addsd	%xmm0, %xmm2
	movabsq	$4835703278458516699, %rcx
	movq	%rdx, -168(%rbp)
	mulq	%rcx
	shrq	$18, %rdx
	cvtsi2sdq	%rdx, %xmm0
	mulsd	%xmm2, %xmm0
	movsd	-128(%rbp), %xmm2
	movsd	-120(%rbp), %xmm3
	subsd	%xmm3, %xmm2
	divsd	%xmm2, %xmm0
	movq	___stderrp@GOTPCREL(%rip), %rdx
	movq	(%rdx), %rdi
	movb	$2, %r8b
	leaq	L_.str1(%rip), %r9
	movq	%rsi, -176(%rbp)
	movq	%r9, %rsi
	movsd	%xmm0, -184(%rbp)
	movapd	%xmm2, %xmm0
	movsd	-184(%rbp), %xmm2
	movsd	%xmm1, -192(%rbp)
	movapd	%xmm2, %xmm1
	movb	%r8b, %al
	movq	%rcx, -200(%rbp)
	movq	%rdx, -208(%rbp)
	callq	_fprintf
	movsd	-128(%rbp), %xmm0
	movsd	-120(%rbp), %xmm1
	subsd	%xmm1, %xmm0
	movq	-72(%rbp), %rax
	movq	%rax, %rcx
	shrq	$32, %rcx
	movq	-176(%rbp), %rdx
	leaq	(%rcx,%rdx), %rcx
	movd	%rcx, %xmm1
	movsd	-192(%rbp), %xmm2
	subsd	%xmm2, %xmm1
	movl	%eax, %eax
	movq	-168(%rbp), %rcx
	leaq	(%rax,%rcx), %rax
	movd	%rax, %xmm3
	addsd	%xmm1, %xmm3
	mulsd	%xmm3, %xmm0
	movq	-64(%rbp), %rax
	movq	-56(%rbp), %rsi
	movq	-200(%rbp), %rdi
	mulq	%rdi
	shrq	$18, %rdx
	imulq	%rsi, %rdx
	movq	%rdx, %rax
	shrq	$32, %rax
	movq	-176(%rbp), %rsi
	leaq	(%rax,%rsi), %rax
	movd	%rax, %xmm1
	subsd	%xmm2, %xmm1
	movl	%edx, %eax
	movl	%eax, %eax
	leaq	(%rax,%rcx), %rax
	movd	%rax, %xmm3
	addsd	%xmm1, %xmm3
	divsd	%xmm3, %xmm0
	movsd	LCPI3_2(%rip), %xmm1
	movapd	%xmm0, %xmm3
	subsd	%xmm1, %xmm3
	cvttsd2siq	%xmm3, %rax
	movabsq	$-9223372036854775808, %rdx
	xorq	%rdx, %rax
	ucomisd	%xmm1, %xmm0
	cvttsd2siq	%xmm0, %rdx
	cmovbq	%rdx, %rax
	movq	%rax, -88(%rbp)
	movsd	-128(%rbp), %xmm0
	movsd	-120(%rbp), %xmm1
	subsd	%xmm1, %xmm0
	movq	-72(%rbp), %rax
	movq	%rax, %rdx
	shrq	$32, %rdx
	leaq	(%rdx,%rsi), %rdx
	movd	%rdx, %xmm1
	subsd	%xmm2, %xmm1
	movl	%eax, %eax
	leaq	(%rax,%rcx), %rax
	movd	%rax, %xmm3
	addsd	%xmm1, %xmm3
	mulsd	%xmm3, %xmm0
	movq	-64(%rbp), %rax
	movq	-56(%rbp), %rdx
	movq	%rdx, %r8
	shrq	$32, %r8
	leaq	(%r8,%rsi), %rsi
	movd	%rsi, %xmm1
	subsd	%xmm2, %xmm1
	movl	%edx, %edx
	leaq	(%rdx,%rcx), %rcx
	movd	%rcx, %xmm2
	addsd	%xmm1, %xmm2
	mulq	%rdi
	movq	%rdx, %rcx
	shrq	$18, %rcx
	cvtsi2sdq	%rcx, %xmm1
	mulsd	%xmm2, %xmm1
	divsd	%xmm1, %xmm0
	movsd	%xmm0, -144(%rbp)
	movq	-208(%rbp), %rcx
	movq	(%rcx), %rcx
	movl	-88(%rbp), %esi
	movb	$1, %dil
	leaq	L_.str2(%rip), %r8
	movb	%dil, -209(%rbp)
	movq	%rcx, %rdi
	movl	%esi, -216(%rbp)
	movq	%r8, %rsi
	movl	-216(%rbp), %ecx
	movl	%ecx, %edx
	movb	-209(%rbp), %cl
	movb	%cl, %al
	callq	_fprintf
	movl	$0, -24(%rbp)
	movl	$0, -20(%rbp)
	movl	-20(%rbp), %eax
	addq	$224, %rsp
	popq	%rbp
	ret
Leh_func_end3:

	.section	__TEXT,__literal8,8byte_literals
	.align	3
LCPI4_0:
	.quad	4517329193108106637
	.section	__TEXT,__text,regular,pure_instructions
	.globl	_wallclock_time
	.align	4, 0x90
_wallclock_time:
Leh_func_begin4:
	pushq	%rbp
Ltmp8:
	movq	%rsp, %rbp
Ltmp9:
	subq	$48, %rsp
Ltmp10:
	leaq	-32(%rbp), %rax
	movabsq	$0, %rcx
	movq	%rax, %rdi
	movq	%rcx, %rsi
	callq	_gettimeofday
	movl	_base.3629(%rip), %eax
	cmpl	$0, %eax
	jne	LBB4_2
	movq	-32(%rbp), %rax
	movq	%rax, _b_val.3627(%rip)
	movl	-24(%rbp), %eax
	movl	%eax, _b_val.3627+8(%rip)
	movl	$1, _base.3629(%rip)
	movabsq	$0, %rax
	cvtsi2sdq	%rax, %xmm0
	movsd	%xmm0, -16(%rbp)
	jmp	LBB4_3
LBB4_2:
	movq	-32(%rbp), %rax
	movq	_b_val.3627(%rip), %rcx
	subq	%rcx, %rax
	cvtsi2sdq	%rax, %xmm0
	movl	-24(%rbp), %eax
	cvtsi2sd	%eax, %xmm1
	movl	_b_val.3627+8(%rip), %eax
	cvtsi2sd	%eax, %xmm2
	subsd	%xmm2, %xmm1
	movsd	LCPI4_0(%rip), %xmm2
	mulsd	%xmm2, %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -40(%rbp)
	movsd	-40(%rbp), %xmm0
	movsd	%xmm0, -16(%rbp)
LBB4_3:
	movsd	-16(%rbp), %xmm0
	movsd	%xmm0, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	addq	$48, %rsp
	popq	%rbp
	ret
Leh_func_end4:

	.globl	_wallclock_time_
	.align	4, 0x90
_wallclock_time_:
Leh_func_begin5:
	pushq	%rbp
Ltmp11:
	movq	%rsp, %rbp
Ltmp12:
	subq	$16, %rsp
Ltmp13:
	callq	_wallclock_time
	movsd	%xmm0, -16(%rbp)
	movsd	-16(%rbp), %xmm0
	movsd	%xmm0, -8(%rbp)
	movsd	-8(%rbp), %xmm0
	addq	$16, %rsp
	popq	%rbp
	ret
Leh_func_end5:

	.section	__TEXT,__cstring,cstring_literals
L_.str:
	.asciz	 "ptr =%d\n"

	.align	3
L_.str1:
	.asciz	 "float pow time=%f gives %f Mflop/s\n"

	.align	3
L_.str2:
	.asciz	 "float pow used %d cycles per instruction %f\n"

.zerofill __DATA,__bss,_base.3629,4,2
.zerofill __DATA,__bss,_b_val.3627,16,3
	.section	__TEXT,__eh_frame,coalesced,no_toc+strip_static_syms+live_support
EH_frame0:
Lsection_eh_frame:
Leh_frame_common:
Lset0 = Leh_frame_common_end-Leh_frame_common_begin
	.long	Lset0
Leh_frame_common_begin:
	.long	0
	.byte	1
	.asciz	 "zR"
	.byte	1
	.byte	120
	.byte	16
	.byte	1
	.byte	16
	.byte	12
	.byte	7
	.byte	8
	.byte	144
	.byte	1
	.align	3
Leh_frame_common_end:
	.globl	_TAB.eh
_TAB.eh:
Lset1 = Leh_frame_end1-Leh_frame_begin1
	.long	Lset1
Leh_frame_begin1:
Lset2 = Leh_frame_begin1-Leh_frame_common
	.long	Lset2
Ltmp14:
	.quad	Leh_func_begin1-Ltmp14
Lset3 = Leh_func_end1-Leh_func_begin1
	.quad	Lset3
	.byte	0
	.byte	4
Lset4 = Ltmp0-Leh_func_begin1
	.long	Lset4
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset5 = Ltmp1-Ltmp0
	.long	Lset5
	.byte	13
	.byte	6
	.align	3
Leh_frame_end1:

	.globl	_Victim.eh
_Victim.eh:
Lset6 = Leh_frame_end2-Leh_frame_begin2
	.long	Lset6
Leh_frame_begin2:
Lset7 = Leh_frame_begin2-Leh_frame_common
	.long	Lset7
Ltmp15:
	.quad	Leh_func_begin2-Ltmp15
Lset8 = Leh_func_end2-Leh_func_begin2
	.quad	Lset8
	.byte	0
	.byte	4
Lset9 = Ltmp2-Leh_func_begin2
	.long	Lset9
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset10 = Ltmp3-Ltmp2
	.long	Lset10
	.byte	13
	.byte	6
	.align	3
Leh_frame_end2:

	.globl	_main.eh
_main.eh:
Lset11 = Leh_frame_end3-Leh_frame_begin3
	.long	Lset11
Leh_frame_begin3:
Lset12 = Leh_frame_begin3-Leh_frame_common
	.long	Lset12
Ltmp16:
	.quad	Leh_func_begin3-Ltmp16
Lset13 = Leh_func_end3-Leh_func_begin3
	.quad	Lset13
	.byte	0
	.byte	4
Lset14 = Ltmp5-Leh_func_begin3
	.long	Lset14
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset15 = Ltmp6-Ltmp5
	.long	Lset15
	.byte	13
	.byte	6
	.align	3
Leh_frame_end3:

	.globl	_wallclock_time.eh
_wallclock_time.eh:
Lset16 = Leh_frame_end4-Leh_frame_begin4
	.long	Lset16
Leh_frame_begin4:
Lset17 = Leh_frame_begin4-Leh_frame_common
	.long	Lset17
Ltmp17:
	.quad	Leh_func_begin4-Ltmp17
Lset18 = Leh_func_end4-Leh_func_begin4
	.quad	Lset18
	.byte	0
	.byte	4
Lset19 = Ltmp8-Leh_func_begin4
	.long	Lset19
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset20 = Ltmp9-Ltmp8
	.long	Lset20
	.byte	13
	.byte	6
	.align	3
Leh_frame_end4:

	.globl	_wallclock_time_.eh
_wallclock_time_.eh:
Lset21 = Leh_frame_end5-Leh_frame_begin5
	.long	Lset21
Leh_frame_begin5:
Lset22 = Leh_frame_begin5-Leh_frame_common
	.long	Lset22
Ltmp18:
	.quad	Leh_func_begin5-Ltmp18
Lset23 = Leh_func_end5-Leh_func_begin5
	.quad	Lset23
	.byte	0
	.byte	4
Lset24 = Ltmp11-Leh_func_begin5
	.long	Lset24
	.byte	14
	.byte	16
	.byte	134
	.byte	2
	.byte	4
Lset25 = Ltmp12-Ltmp11
	.long	Lset25
	.byte	13
	.byte	6
	.align	3
Leh_frame_end5:


.subsections_via_symbols
