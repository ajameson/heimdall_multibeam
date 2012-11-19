
.cu.o:
	$(CUDA_NVCC) $(CUDA_CFLAGS) $(CUDA_NVCC_CFLAGS) -o $@ -c $<

.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(CUDA_NVCC) $(CUDA_CFLAGS) $(CUDA_NVCC_CFLAGS) --compiler-options=\" $(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" -c $<
