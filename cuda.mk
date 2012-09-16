#$(CUDA_NVCC) -gencode=arch=compute_20,code=sm_20 -o $@ -c $<
#$(top_srcdir)/cudalt.py $@ $(CUDA_NVCC) -gencode=arch=compute_20,code=sm_20 --compiler-options=\" $(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" -c $<

.cu.o:
	$(CUDA_NVCC) -O3 -arch sm_20 -o $@ -c $<

.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(CUDA_NVCC) -O3 -arch sm_20 --compiler-options=\" $(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" -c $<
