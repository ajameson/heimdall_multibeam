/***************************************************************************
 *
 *   Copyright (C) 2012 by Ben Barsdell and Andrew Jameson
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/
:wn

cpp_objects: $(CPP_SOURCES) $(HEADERS)
	$(GXX) -c $(GXX_FLAGS) $(INCLUDE) $(CPP_SOURCES)

cuda_objects: $(CUDA_SOURCES) $(HEADERS)
	$(NVCC) -c $(NVCC_FLAGS) $(INCLUDE) $(CUDA_SOURCES)

// 8-bit integer out-of-place version
void repartition(const unsigned char* in, size_t count,
                 unsigned char* out) {
	MPI_Comm comm = MPI_COMM_WORLD;
	int proc_count;
	MPI_Comm_size(comm, &proc_count);
	int proc_idx;
	MPI_Comm_rank(comm, &proc_idx);
	
	std::vector<int> send_counts(proc_count);
	std::vector<int> send_offsets(proc_count);
	std::vector<int> recv_counts(proc_count);
	std::vector<int> recv_offsets(proc_count);
	send_offsets[0] = 0;
	recv_offsets[0] = 0;
	for( size_t i=1; i<proc_count; ++i ) {
		send_counts[i] = count / proc_count + (i < in_count % proc_count);
		recv_counts[i] = send_counts[i];
	}
	for( size_t i=1; i<proc_count; ++i ) {
		send_offsets[i] = send_offsets[i-1] + send_counts[i-1];
		recv_offsets[i] = send_offsets[i];
	}
	
	MPI_Alltoallv((void*)in, &send_counts[0], &send_offsets[0], MPI_BYTE,
	              (void*)out, &recv_counts[0], &recv_offsets[0], MPI_BYTE,
	              comm);
}


// 8-bit integer version
void repartition(unsigned char* data, size_t count) {
	MPI_Comm comm = MPI_COMM_WORLD;
	int proc_count;
	MPI_Comm_size(comm, &proc_count);
	int proc_idx;
	MPI_Comm_rank(comm, &proc_idx);
	
	std::vector<int> recv_counts(proc_count);
	std::vector<int> recv_offsets(proc_count);
	recv_offsets[0] = 0;
	for( size_t i=1; i<proc_count; ++i ) {
		recv_counts[i] = count / proc_count + (i < in_count % proc_count);
	}
	for( size_t i=1; i<proc_count; ++i ) {
		recv_offsets[i] = recv_offsets[i-1] + recv_counts[i-1];
	}
	
	MPI_Alltoallv((void*)MPI_IN_PLACE, 0, 0, 0,
	              (void*)&data[0], &recv_counts[0], &recv_offsets[0], MPI_BYTE,
	              comm);
}


// Explicit template instantiations
// TODO: This isn't ideal, but otherwise we run into trouble separating
//         the normal/CUDA/MPI sources for their respective compilers.
//template void repartition<unsigned char>(unsigned char* data, size_t count);
//template void repartition<float>(float* data, size_t count);

#repartition.o: repartition.cpp $(HEADERS)
#	$(GXX) -c $(MPICXX_FLAGS) $(INCLUDE) repartition.cpp

/*
		  If you divide dm_count by nbeams and round up and pad, then the last
		    process may end up with padding at intervals through its dst array.
		  Actually, rounding up is not the best strategy due to balancing issues.
		    An optimal partitioning is n[i] = N/P + i<N%P
		  
		  
		*/

/*
		char name[MPI_MAX_PROCESSOR_NAME];
		int len;
		MPI_Get_processor_name(name, &len);
		cout << "I'm running on " << name << endl;
		*/

//cout << "dm = " << dm << endl;
				//cout << "local_dm_count = " << local_dm_count << endl;
				//cout << "nsamps_computed = " << nsamps_computed << endl;
				//cout << "offset = " << offset << endl;

/*
	MPI_Alltoallv((void*)MPI_IN_PLACE, 0, 0, 0,
	              (void*)&data[0],
	              &recv_counts[0], &recv_offsets[0],
	              get_mpi_datatype(T()),
	              comm);
	*/
	// WAR for lack of MPI_IN_PLACE support in OpenMPI's Alltoall (2012-26-01)

		//cout << "("<<proc_idx<<")"<<"send_counts[" << i << "] = " << send_counts[i] << endl;
		//cout << "("<<proc_idx<<")"<<"recv_counts[" << i << "] = " << recv_counts[i] << endl;
		//cout << "("<<proc_idx<<")"<<"send_offsets[" << i << "] = " << send_offsets[i] << endl;
		//cout << "("<<proc_idx<<")"<<"recv_offsets[" << i << "] = " << recv_offsets[i] << endl;


	thrust::exclusive_scan(d_in_begin, d_in_end, d_scanned.begin()+w/2);
	// Now we need to fill the ends with 
	thrust::fill(d_scanned.begin(), d_scanned.begin()+w/2, hd_float(0));
	thrust::fill(d_scanned.end()-w/2, d_scanned.end, d_scanned[count-w/2]);
	
/*
	// TODO: Scan by key and then use lower_bound to find splitting elements
	//         Can interpolate between with (x-x0)/(x1-x0)
	thrust::device_vector<unsigned int> d_scanned_giant_inds(d_giant_inds.size());
	thrust::inclusive_scan_by_key(d_giant_inds.begin(),
	                              d_giant_inds.end(),
	                              d_giant_inds.begin(),
	                              d_scanned_giant_inds.begin(),
	                              adjacency<unsigned int>());
	thrust::lower_bound_by_key(
	*/
	
		              //d_keys_result.begin(), // WAR for not having discard_iterator

	//thrust::device_vector<hd_size>  d_keys_result(d_giant_data.size());

/*
	thrust::transform_if(make_counting_iterator<unsigned int>(0),
	                     make_counting_iterator<unsigned int>(count),
	                     d_data_begin, // the stencil
	                     d_data_begin,
	                     clip<hd_float>(clip_limit),
	*/

/*
	thrust::transform(d_data_begin, d_data_end,
	                  d_scanned.begin()+r,
	                  clip<hd_float>(clip_limit));
	*/


				//hd_float* in = thrust::raw_pointer_cast(&d_beam_series[0]);
				//hd_float* out = thrust::raw_pointer_cast(&d_filtered_series[0]);
				

2-bit
unpack/extend to at least 8-bit
scrunch_x2

2-bit (3/3)
scrunch --> 4-bit (6/15)
scrunch (12/15)
scrunch --> 8-bit (24/256)
scrunch (48/256)
scrunch (96/256)
scrunch (192/256)
scrunch --> 16-bit (384/65536)

OR

2-bit (3/3)
scrunch --> 8-bit (6/256)
scrunch (12/256)
scrunch (24/256)
scrunch (48/256)
scrunch (96/256)
scrunch (192/256)


		/*
		  if( time to scrunch ) {
		      h_scrunched.resize(nsamps_computed/2 * pl->params.nchans);
		      scrunch_filterbank_x2(h_filterbank_expanded --> h_scrunched_filterbank,
		                            in_nbits);
		      swap(h_filterbank_expanded, h_scrunched_filterbank);
		      nsamps_computed /= 2;
		      dt *= 2;
		      
		  }
		  dedisperse h_filterbank_expanded --> h_dm_series;
		*/
		
		hd_float dm0      = dm_list[gulp_dm];
		hd_float delta_dm = dm_list[gulp_dm+1] - dm0;
		hd_float smearing = get_smearing(pl->params.dt,
		                                 pl->params.dm_pulse_width*1e-6,
		                                 pl->params.f0,
		                                 pl->params.nchans, pl->params.df,
		                                 dm0, delta_dm);
		// The smearing after time scrunching
		hd_float smearing2 = get_smearing(2 * pl->params.dt,
		                                  pl->params.dm_pulse_width*1e-6,
		                                  pl->params.f0,
		                                  pl->params.nchans, pl->params.df,
		                                  dm0, delta_dm);
		
		if( smearing2 / smearing < pl->params.scrunch_tol ) {
			start_timer(filter_timer);
			scrunch_filterbank_x2(h_filterbank, nsamps, nbits, h_filterbank);
			nsamps_computed /= 2;
			pl->params.dt *= 2;
			stop_timer(filter_timer);
		}
		
		/*
		if( rank == 0 ) {
			cout << "*** smearing/dt        = " << smearing / pl->params.dt << endl;
			cout << "    smearing2/smearing = " << smearing2/smearing << endl;
		}
		*/
		//if( smearing > 2 * pl->params.dt ) {
			
		//}
		

// First we expand to a minimum number of bits per sample, so that the
	//   filterbank can be scrunched at will without major loss of precision.
	hd_size min_nbits = 8;
	if( nbits < min_nbits ) {
		hd_size size_bytes = nsamps * pl->params.nchans * nbits / 8;
		//h_filterbank_tmp.resize(size_bytes * (min_nbits / nbits));
		expand_nbits(h_filterbank, size_bytes, nbits,
		             //&h_filterbank_tmp[0],
		             h_filterbank,
		             min_nbits);
		//h_filterbank = &h_filterbank_tmp[0];
		nbits = min_nbits;
	}
	


	//host_vector<hd_byte> h_filterbank_tmp;
	//host_vector<hd_byte> h_scrunched_filterbank();
	

	/*
private:
	
	// TODO: NEED PROPER COPY/ASSIGN OPS DUE TO DYN MEM ALLOC!
	//       See http://drdobbs.com/cpp/205918714?pgno=3 for a really nice
	//         solution.
	
	MatchedFilterPlan_impl* m_impl;
	*/

		//dedisp_set_dm_list(pl->dedispersion_plan, &dm_list[gulp_dm], dm_gulp_size);

// TODO: Don't think it's necessary. Just forgot to increase dt.
				//nsamps_smooth = min(nsamps_smooth, cur_nsamps);
//	thrust::raw_pointer_cast(&d_beam_series[beam *
				//	                                        series_stride]);

			// We must remember the actual cur_nsamps here before scrunching
			//hd_size unscrunched_cur_nsamps = cur_nsamps;
			
/*
				// Swap the pointers for double-buffering
				std::swap(beam_series, filtered_series);
				
				start_timer(filter_timer);
				if( filter_width != 1 ) {
					scrunch_x2(beam_series, cur_nsamps*beam_count,
					           filtered_series);
					cur_nsamps /= 2;
				}
				stop_timer(filter_timer);
				*/

//thrust::raw_pointer_cast(&d_beam_series[beam *
					//                                        series_stride]);
	//thrust::raw_pointer_cast(&d_filtered_series[beam *
						//                                            series_stride]);
					hd_float* beam_ptr = &beam_series[beam*series_stride];
					// Note: This modifies cur_nsamps (/2)
					//error = scrunch_x2(in, cur_nsamps, filtered);


					/*
					  start_timer(normalise_timer);
					  error = normalise(filtered, cur_nsamps);
					  stop_timer(normalise_timer);
					  // TESTING
					  //error = normalise2(filtered, cur_nsamps);
					  if( error != HD_NO_ERROR ) {
					  return throw_error(error);
					  }
					*/

//	thrust::raw_pointer_cast(&d_filtered_series[beam *
					//	                                            cur_nsamps]);



				// Swap the pointers for double-buffering
				//std::swap(beam_series, filtered_series);
				
			//cur_nsamps = unscrunched_cur_nsamps;

/*			
Apply time-scrunching at certain DM (gulp) intervals and then
        go back to using proper convolutions for the filtering.
        This will save an overall factor of ~2x while maximising
          sensitivity.

We can directly swap( ) vectors rather than swapping pointers to
        them if we want to. Check if this is worthwhile.
      
IMPORTANT: Within dedisp, we really want to generate a DM list AND a dt list!
             Could toggle time-scrunching functionality on/off
             Could provide a tolerance factor for generating the dt list
             Would have to break up the computed DMs into blocks of equal dt
             The only major problem is scrunching when nbits < 8
               Could expand into a tmp array...
               ACTUALLY, from a quick Python test, it looks like it doesn't really
                 matter! And really, if the user wants to they can just expand
                 it themselves before giving it to dedisp!

 11100100 3210
+00011011 0123
=11111111 3333 /2 =1111
/2
=01111111

*/


// TODO: Where should this go?
hd_float get_smearing(hd_float dt, hd_float pulse_width,
                      hd_float f0, hd_size nchans, hd_float df,
                      hd_float DM, hd_float deltaDM) {
	hd_float W         = pulse_width;
	hd_float BW        = nchans * abs(df);
	hd_float fc        = f0 - BW/2;
	hd_float inv_fc3   = 1./(fc*fc*fc);
	hd_float t_DM      = 8.3*BW*DM*inv_fc3;
	hd_float t_deltaDM = 8.3/4*BW*nchans*deltaDM*inv_fc3;
	hd_float t_smear   = sqrt(dt*dt + W*W + t_DM*t_DM + t_deltaDM*t_deltaDM);
	return t_smear;
}

// Read and apply zero dm RFI mask
	//read_zero_dm_rfi_mask(&zero_dm_rfi_mask[0], ...);
	//zap_zero_dm_rfi(&filterbank[0], &zero_dm_rfi_mask[0], ...);
	
	/*
	  Dedisperse @ DM=0 USING CUSTOM METHOD (i.e., not dedisp)
	  // MPI_Gather the time series on proc 0
	  MPI_Gather(&dm0_series[0], cur_nsamps, get_mpi_datatype(OutType),
	             &dm0_beam_series[0], cur_nsamps, get_mpi_datatype(OutType),
	             0, MPI_COMM_WORLD);
	  If proc 0:
          Execute create_multibeam_coinc_mask(&dm0_beam_series[0], ..., &mask[0], ...)
	  MPI_Scatter the mask to all the other procs
	  MPI_Scatter(&mask[0], cur_nsamps, get_mpi_datatype(OutType),
	              &mask[0], cur_nsamps, get_mpi_datatype(OutType),
	              0, MPI_COMM_WORLD);
	  Execute zap_filterbank_rfi(&mask[0], ...)
	*/
	//dedisp_error derror;
	// --------------



// From ccl_test.py
// ----------------

    """
    new_coords, old_labels, new_labels = find_adjacency(new_coords, \
                                                            old_labels, new_labels, \
                                                            x_fastest, adjacent_x)
    #print "1)New coords:"
    #print new_coords
    #print "1)New labels:"
    #print new_labels
    #print "1)Corresponding old labels:"
    #print old_labels
    
    new_coords, old_labels, new_labels = find_adjacency(new_coords, \
                                                            old_labels, new_labels, \
                                                            y_fastest, adjacent_y)
    
    print "2)New coords:"
    print new_coords
    print "2)New labels:"
    print new_labels
    print "2)Corresponding old labels:"
    print old_labels
    """
    

# TESTING
#a = [1, 2, 3]
#a[1:] = [4,9]
#print a
#a[1:] = map(lambda x:x*x, a[1:])
#print a

def adjacent_x(b, a):
    return b[1] == a[1] and b[0] == a[0] + 1
def adjacent_y(b, a):
    return b[0] == a[0] and b[1] == a[1] + 1

def find_adjacency(coords, old_labels, labels, key_, adjacency):
    coords_olds_labels = zip(coords, old_labels, labels)
    #print "coords_labels:"
    #print coords_labels
    
    # TODO: Make this function return the results in the same order as the
    #         input so that we don't have to worry about tracking the old_labels.
    
    sorted_coords_olds_labels = sorted(coords_olds_labels, \
                                           key=lambda (coord,old_label,label): \
                                           key_(coord))
    coords, old_labels, labels = zip(*sorted_coords_olds_labels)
    #print "Sorted coords_olds_labels:"
    #print sorted_coords_olds_labels
    
    sorted_coords_labels = zip(coords, labels)
    #new_coords_labels = [sorted_coords_labels[0]]
    new_coords_labels = sorted_coords_labels
    
    new_coords_labels[1:] = map(lambda b,a: b if not adjacency(b[0], a[0]) \
                                    else (b[0],min(a[1],b[1])), \
                                    new_coords_labels[1:], \
                                    new_coords_labels[:-1])
    #new_coords_labels[1:] = new2_coords_labels
    #new_coords, new_labels = zip(*new_coords_labels)
    #print "  1)Coords, old and new labels:"
    #print " ", coords
    #print " ", old_labels
    #print " ", new_labels
    new_coords_labels[:-1] = map(lambda a,b: a if not adjacency(b[0], a[0]) \
                                 else (a[0],min(a[1],b[1])), \
                                 new_coords_labels[:-1], \
                                 new_coords_labels[1:])
    #new_coords_labels[:-1] = new2_coords_labels
    new_coords, new_labels = zip(*new_coords_labels)
    #print "  2)Coords, old and new labels:"
    #print " ", coords
    #print " ", old_labels
    #print " ", new_labels
    return (new_coords, old_labels, new_labels)

def exclusive_scan(values, init):
    result = []
    sum = init
    for x in values:
        result.append(sum)
        sum += x
    return result


"""
it_count = 0
for i in range(1,it_count):
    coords, labels = iterate_ccl(coords, labels)
    print "Iteration %i:" % (i+1)
    print zip(coords, labels)
    print labels
"""


// ----------------


	/*
	t mask data exscan 
	0   0  ABCD   0     
	1   1  EFGH   0     
	2   0  IJKL   1     
	3   0  MNOP   1     
	4   1  QRST   1     
	5   1  UVWX   2     
	  
	  t = range(mask);
	  //(good, bad) = separate(t, mask==0);
	  good = t(mask==0);
	  bad  = t(mask==1);
	  sub  = shuffle(good)(0:len(bad));
	  c = range(nchans);
	  out(c,bad) = in(c,sub);
	  
	  j = range(bad);
	  out[c + bad[j] * nchans] = in[c + sub[j] * nchans];
	  c = i % nchans;
	  j = i / nchans;
	  
	*/
	/*
	make_permutation_iterator(make_transform_iterator(make_counting_iterator(0),
	                                                  _1 % val(nchans) +
	                                                  _1 / val(nchans) 
	                                                  
	If we did this inside dedisp instead then we could do it after transposition
	TODO: The below code outline should work, but still requires dealing with nbits.
	        Could perhaps just expand it to floats as we do each channel.
	for each channel c:
	  get_baseline(channel c time series, means);
	  rms = get_rms(channel c time series);
	  transform_if(make_zip_iterator(make_tuple(make_counting_iterator(0),
	                                            means)),
	               make_zip_iterator(make_tuple(make_counting_iterator(0),
	                                            means)) + mask_nsamps,
	               mask,
	               channel c time series begin,
	               AddGaussianNoise(rms),
	               thrust::identity<int>());
	                                                  
	*/
	
	/*
	  
	  in:     01 10 11 02 12 20 21
	  sorted: 01 02 10 11 12 20 21
	          
	*/
	
// TODO: Zap the whole band for each masked time sample
	//         Easiest way to do this is to grab random nearby samples
	//           from the same band.
	/*
	  
	  for each masked time sample t:
	      for each channel c:
	          hd_size t_min = max(t-smooth_radius, 0);
	          hd_size t_max = min(t+smooth_radius, nsamps);
	          rand_t = t_min + rand() * (t_max-t_min);
              out[t,c] = in[rand_t,c];
	  
	 */

// From matched_filter.cu
// ----------------------
//thrust::exclusive_scan(m_scanned.begin(), m_scanned.end(),
		//                       m_scanned.begin());
/*
		thrust::transform(m_scanned.begin()+width, m_scanned.end(),
		                  m_scanned.begin(),
		                  d_out_begin,
		                  averager<hd_float>(width));
		*/
// ----------------------

// TODO: Remove this when happy with the OOP version below
					//error = matched_filter(beam_ptr, cur_nsamps,
					//                       filter_width, filtered);
					// Note: When filter_width==1 we could just do a simple copy,
					//         but it's probably no problem to just leave it.
					// HACK TESTING
					//hd_float* beam_ptr = &beam_series[beam*series_stride];
					/*
					if( filter_width == 1 ) {
						thrust::copy(thrust::device_ptr<hd_float>(beam_ptr),
						             thrust::device_ptr<hd_float>(beam_ptr) + cur_nsamps,
						             thrust::device_ptr<hd_float>(filtered));
					}
					else {
					*/

/*
			// TODO: Merge this loop with the previous one
			// Normalise each beam
			for( hd_size beam=0; beam<beam_count; ++beam ) {
				hd_float* beam_ptr = &beam_series[beam*series_stride];
				
			}
			
			*/
			

// From matched_filter.h
// ---------------------

// TODO: Remove this when satisfied with the OOP version below
hd_error matched_filter(const hd_float* d_in,
                        hd_size         count,
                        hd_size         filter_width,
                        hd_float*       d_out);

struct MatchedFilterPlan_impl;

struct MatchedFilterPlan {
	MatchedFilterPlan();
	hd_error prep(const hd_float* d_in, hd_size count, hd_size max_width);
	hd_error exec(hd_size width, hd_float* d_out);
	
private:
	boost::shared_ptr<MatchedFilterPlan_impl> m_impl;
};

// ---------------------

/*
	std::cout << "labels:" << std::endl;
	thrust::copy(d_labels_begin, d_labels_begin+count,
	             std::ostream_iterator<hd_size>(std::cout, "\t"));
	std::cout << std::endl;
	*/
	/*
	labels = sequence(count);
	search_offsets = [(0,-1,0), (0,-1,-1), (0,0,-1)];
	for offset in search_offsets:
	    search_vals = project(begins + (-r+1, offset));
	    neib_begins = lower_bound(project(ends), search_vals);
	    search_vals = project(begins + (+r, offset));
	    neib_ends   = upper_bound(project(ends), search_vals);
	    
	    new_labels = transform(labels,
	                           zip(neib_begins, neib_ends),
	                           range_min(labels));
	    labels = new_labels;
	
	for i in range(count):
        j = new_labels[i]
        while new_labels[j] != j:
            j = new_labels[j]
        new_labels[i] = j
	*/
	
/*
			  events[filter];
					  
			  0123 567
			  0 2   6 8
			  00122356678 merged
			  0 2   6     flagged
					  
			  for( hd_size i=0; i<filter_count; ++i ) {
					      
			  }
					  
			  e e e e   e e e
			  e   e       e   e
			  e e   e e e e   e e
					  
			*/

//if( new_cand_count > 0 ) {
			//cout << "FOUND " << new_cand_count << " NEW CANDIDATES" << endl;
			//}
			//if( new_cand_count > 0 ) {
			/*for( hd_size nearby_beam=0; nearby_beam<3; ++nearby_beam ) {
			  std::stringstream ss;
			  hd_size b = (beam + nearby_beam) % beam_count;
			  ss << "candidate"
			  << "_snr" << h_giant_peaks[prev_giant_count]
			  << "_beam" << b+1
			  << "_dm" << actual_dm
			  << "_filter" << filter_width
			  << "_samp" << h_giant_inds[prev_giant_count]
			  << ".tim";
			  //cerr << "Writing candidate to " << ss.str() << endl;
			  //hd_float* p =&filtered_series[b*cur_nsamps];
			  //	thrust::raw_pointer_cast(&d_filtered_series[b *
			  //                                            cur_nsamps]);
							
			  write_device_time_series(p,
			  cur_nsamps,
			  cur_dt,
			  ss.str().c_str());
							
			  }*/
			//cout << "PROC_COUNT = " << proc_count << endl;
					
			//if( rank == 0 ) {
			// HACK to emulate an atomic write
			std::ofstream cand_file("candidates.dat", std::ios::app);
			while( !cand_file.is_open() ) {
				cand_file.open("candidates.dat", std::ios::app);
			}
			for( hd_size i=prev_giant_count; i<h_giant_peaks.size(); ++i ) {
				// HACK %13
				cand_file << h_giant_peaks[i]
					// TODO: This sample index may need adjusting for the filtered offset
				          << "\t" << first_idx+h_giant_inds[i]*cur_scrunch+cur_filtered_offset
				          << "\t" << (h_giant_beams[i] + pl->params.first_beam)%13 + 1
				          << "\t" << h_giant_filters[i]
				          << "\t" << h_giant_dms[i]
				          << "\t" << h_giant_members[i] << "\n";
			}
			cand_file.close();
			//}
					
			//for (int p=0; p<proc_count; p++) {
			//if( p == rank ) {
			//cout << "WRITING CANDIDATES FROM PROCESS " << p << endl;
								
			//}
			//MPI_Barrier(MPI_COMM_WORLD);
			//}
			/*
			  if( rank == 0 ) {
			  std::ofstream cand_file("candidates.dat", std::ios::app);
			  for( hd_size i=prev_giant_count; i<h_giant_peaks.size(); ++i ) {
			  cand_file << h_giant_peaks[i]
			  << "\t" << h_giant_inds[i]
			  << "\t" << h_giant_beams[i]
			  << "\t" << h_giant_filters[i]
			  << "\t" << h_giant_dms[i] << "\n";
			  }
			  cand_file.close();
			  }
			*/
			//}
					
			////report_candidates(&candidates[0]);
				


	// TESTING
	/*
	  if( pl->params.verbosity >= 3 ) {
	  for( hd_size i=0; i<h_giant_peaks.size(); ++i ) {
	  cout << "Found candidate with SNR " << h_giant_peaks[i]
	  << "\tat index " << first_idx+h_giant_inds[i]*scrunch_factors[???]
	  << "\tin beam " << h_giant_beams[i]+1 + pl->params.first_beam
	  << "\tusing filter width of " << h_giant_filters[i]
	  << "\twith DM of " << h_giant_dms[i]
	  << endl;
	  }
	  }
	*/
	

		
			// TODO: 1) Call find_giants to find all time-connected events
			//            for this beam, DM and filter.
			//       2) Call label_candidate_clusters to produce
			//            'island labels'.
			//       3) Call a big reduce_by_key to find the peak SNR
			//            (and other props) of each island.
			//            The result of this can be output as a candidate list.
			

/*
	thrust::transform(d_group_members.begin(), d_group_members.end(),
	                  make_constant_iterator<hd_size>(min_group_members),
	                  d_group_classes.begin(),
	                  thrust::less_than<hd_size>());
	*/

/*
	thrust::transform_if(d_group_dms.begin(), d_group_dms.end(),
	                     make_constant_iterator<hd_size>(rfi_max_dm),
	                     d_group_classes.begin(),
	                     d_group_classes.begin(),
	                     thrust::less_than<hd_size>(),
	                     thrust::not1(thrust::identity<hd_size>()));
	*/
	

/*
	if( giant_data_count > 0 ) {
		d_giant_data_segments.front() = 1;
	}
	
	// Note: Here we want to explicitly count the first segment; later it
	//         will be implicit.
	hd_size giant_count = thrust::count_if(d_giant_data_segments.begin(),
	                                       d_giant_data_segments.end(),
	                                       thrust::identity<hd_float>());
	*/

/*
	if( giant_count_new != giant_count ) {
		std::cout << "WE HAVE A PROBLEM" << std::endl;
		std::cout << "orig = " << giant_count << std::endl;
		std::cout << "new  = " << giant_count_new << std::endl;
	}
	*/


	/*
	  0001001010001000 head flags
	  0001112233334444 inc scan (wrong)
	  0000111223333444 exc scan (right)
	  
	  0100100101
	  0111222334 inc scan (right)
	  0011122233
	  
	  01001001010
	  01112223344 inc scan (right)
	  00111222334
	  
	  1100100101
	  1222333445
	  0122233344 exc scan (right)
	  
	  11001001010
	  12223334455
	  01222333445 exc scan (right)
	  
	  1001001010001000 head flags
	  1112223344445555 inc scan (wrong)
	  0111222334444555 exc scan (right)
	  
	  scatter_if(data_inds, inc_scan, head_flags, begins);
	 */
	
/*
	// Create an array of the number of members in each giant
	thrust::device_vector<int> d_giant_members(giant_count);
	thrust::transform(d_giant_ends.begin(), d_giant_ends.end(),
	                  d_giant_begins.begin(),
	                  d_giant_members.begin(),
	                  thrust::minus<hd_float>());
	*/
	// Now we compute the FWHM of each giant
	// ---------------------------------------------------
	/*
	// Square the giant data
	thrust::transform(d_giant_data.begin(), d_giant_data.end(),
	                  square<hd_float>());
	// Sum the squares
	reduce_by_key(d_giant_data_inds.begin(), // the keys
	              d_giant_data_inds.end(),
	              d_giant_data.begin(),
	              d_giant_starts.begin(),
	              d_giant_widths.begin(),
	              adjacent<hd_size>(),
	              thrust::plus<hd_float>());
	
	thrust::adjacent_difference(d_giant_starts.begin(),
	                            d_giant_starts.begin()+giant_count,
	                            
	
	thrust::transform(d_giant_widths.begin(), d_giant_widths.end(),
	                  sumsqs_to_stddev(
	*/
	// ---------------------------------------------------
	
	// Append the results to the existing results on the host
	/*
	h_giant_peaks.resize(h_giant_peaks.size()+giant_count);
	h_giant_inds.resize(h_giant_inds.size()+giant_count);
	h_giant_members.resize(h_giant_members.size()+giant_count);
	thrust::copy(d_giant_peaks.begin(), d_giant_peaks.begin()+giant_count,
	             h_giant_peaks.end()-giant_count);
	thrust::copy(d_giant_inds.begin(), d_giant_inds.begin()+giant_count,
	             h_giant_inds.end()-giant_count);
	thrust::copy(d_giant_members.begin(), d_giant_members.begin()+giant_count,
	             h_giant_members.end()-giant_count);
	*/


/*
  
  Sort by x+y*w
  For each pixel x,y:
    tmp = out[x,y]
    if pxl[x+1,y] exists: tmp      = min(pxl[x,y],pxl[x+1,y])
    if pxl[x-1,y] exists: out[x,y] = min(tmp,pxl[x-1,y])
  Sort by y+x*h
  For each pixel x,y:
    tmp = out[x,y]
    if pxl[x,y+1] exists: tmp      = min(pxl[x,y],pxl[x,y+1])
    if pxl[x,y-1] exists: out[x,y] = min(tmp,pxl[x,y-1])
  Sort (in,out) by out
  Flag where in==out
  Fill segments of out with head value (where in==out)
  Sort (in,out) by in
  Repeat starting with out
*/

/*
	  
	  0125689A
	  10010100
	  0  5 8  begins
	    3 7 A
	 */
	

//thrust::device_vector<int> d_giant_begins(giant_count);
	//thrust::device_vector<int> d_giant_ends(giant_count);


	// Note: Here we allocate the max needed instead of pre-computing the size
	//thrust::device_vector<hd_float> d_giant_peaks(d_giant_data.size());
	//thrust::device_vector<hd_size>  d_giant_inds(d_giant_data.size());
	//thrust::device_vector<hd_size>  d_giant_widths(d_giant_data.size());
	//thrust::device_vector<hd_size>  d_giant_starts(d_giant_data.size());
	
/*
	// Finally, trace equivalency chains to find the final labels
	// Note: We do this sequentially on the host as it is more efficient
	thrust::host_vector<hd_size> h_labels(d_labels_begin,
	                                      d_labels_begin+count);
	for( hd_size i=0; i<h_labels.size(); ++i ) {
		hd_size j = h_labels[i];
		while( h_labels[j] != j ) {
			j = h_labels[j];
		}
		h_labels[i] = j;
	}
	// Copy back to the device
	thrust::copy(h_labels.begin(), h_labels.end(),
	             d_labels_begin);
	*/

// Create an array of head flags indicating candidate segments
	thrust::device_vector<int> d_giant_segments(giant_count);
	thrust::adjacent_difference(d_giant_labels.begin(),
	                            d_giant_labels.end(),
	                            d_giant_segments.begin(),
	                            thrust::not_equal_to<hd_size>());
	if( giant_count > 0 ) {
		d_giant_segments[0] = 1;
	}
	hd_size group_count = thrust::count_if(d_giant_segments.begin(),
	                                       d_giant_segments.end(),
	                                       thrust::identity<hd_size>());
	

/*
	thrust::sort_by_key(d_giant_labels.begin(), d_giant_labels.end(),
	                    make_zip_iterator(thrust::make_tuple(d_giant_peaks.begin(),
	                                                         d_giant_inds.begin(),
	                                                         d_giant_begins.begin(),
	                                                         d_giant_ends.begin(),
	                                                         d_giant_filter_inds.begin(),
	                                                         d_giant_dm_inds.begin())));
	*/
	


	hd_size* d_giant_begins_ptr = thrust::raw_pointer_cast(&d_giant_begins[0]);
	hd_size* d_giant_ends_ptr = thrust::raw_pointer_cast(&d_giant_ends[0]);
	hd_size* d_giant_filter_inds_ptr = thrust::raw_pointer_cast(&d_giant_filter_inds[0]);
	hd_size* d_giant_dm_inds_ptr = thrust::raw_pointer_cast(&d_giant_dm_inds[0]);
	hd_size* d_giant_labels_ptr = thrust::raw_pointer_cast(&d_giant_labels[0]);
	


	/*
	label_candidate_clusters(giant_count,
	                         d_giant_begins_ptr,
	                         d_giant_ends_ptr,
	                         time_count,
	                         d_giant_filter_inds_ptr,
	                         filter_count,
	                         d_giant_dm_inds_ptr,
	                         dm_count,
	                         pl->params.cand_sep_time,
	                         d_giant_labels_ptr);
	*/

/*
	                         d_giant_begins_ptr,
	                         d_giant_ends_ptr,
	                         d_giant_filter_inds_ptr,
	                         d_giant_dm_inds_ptr,
	                         */


/*
	int root_proc = 0;
	MPI_Comm comm = MPI_COMM_WORLD;
	std::vector<int> cand_counts;
	std::vector<int> cand_offsets;
	if( rank == root_proc ) {
		cand_counts.resize(beam_count);
		cand_offsets.resize(beam_count);
	}
	MPI_Gather((void*)&giant_count, 1, get_mpi_datatype(giant_count),
	           (void*)&cand_counts[0], 1, get_mpi_datatype(giant_count),
	           root_proc, comm);
	cand_offsets[0] = 0;
	for( hd_size b=1; b<beam_count; ++b ) {
		cand_offsets[b] = cand_offsets[b-1] + cand_counts[b-1];
	}
	hd_size total_cand_count = (cand_offsets[beam_count-1] +
	                            cand_counts[beam_count-1]);
	
	// (1)
	std::vector<hd_float> all_cand_peaks;
	if( rank == root_proc ) {
		all_cand_peaks.resize(total_cand_count);
	}
	MPI_Gatherv((void*)&h_giant_peaks[0], giant_count,
	            get_mpi_datatype(h_giant_peaks[0]),
	            (void*)&all_cand_peaks[0], &cand_counts[0], &cand_offsets[0],
	            get_mpi_datatype(h_giant_peaks[0]),
	            root_proc, comm);
	*/


/*
// TODO: This is just a dummy implementation for testing the pipeline
size_t read_source_data(size_t nsamps, char* data, size_t* nbits) {
	std::ifstream in_file("test.fil", std::ios::binary);
	SigprocHeader header;
	read_header(in_file, header);
	
	*nbits = header.nbits;
	size_t nchan_bytes = header.nchans * header.nbits / (8*sizeof(char));
	
	size_t bytes_read = in_file.readsome((char*)&data[0], nsamps*nchan_bytes);
	in_file.close();
	return bytes_read / nchan_bytes;
}
*/

// From matched_filter.cu
		// We copy the input into the middle of the padded m_scanned array
		// Note: Divide and round up (solves width=odd cases)
		//hd_size offset = (m_max_width-1) / 2 + 1;
		//thrust::copy(d_in_begin, d_in_end, m_scanned.begin()+offset);
		//thrust::copy(d_in_begin, d_in_end, m_scanned.begin());
		// Pre-compute the prefix-sum
		//thrust::inclusive_scan(m_scanned.begin(), m_scanned.end(),
		//                       m_scanned.begin());


			// TODO: Need to account for some overlap between data gulps!
			//         I.e., due to dispersion (and filtering etc.) we need to
			//           re-process some data at the end of the gulp.
			//         This also needs to be considered in the candidate detection.
			//           I.e., we don't want candidates from the overlap region
			
// HACK to emulate an atomic write
	//while( !cand_file.is_open() ) {
	//	cand_file.open("candidates.dat", std::ios::app);
	//}

/*
		cout << "All DM inds:" << endl;
		std::copy(h_all_dm_inds.begin(), h_all_dm_inds.end(),
		          std::ostream_iterator<hd_size>(cout, " "));
		*/

/*
		// Generate a DM histogram from the candidates
		// -------------------------------------------
		hd_float dm_min = pl->params.dm_min;
		hd_float dm_max = pl->params.dm_max;
		hd_size  dm_hist_nbins  = (hd_size)(sqrt(total_cand_count));
		hd_float dm_hist_binsep = (dm_max - dm_min) / dm_hist_nbins;
		std::vector<hd_size> dm_hist(dm_hist_nbins);
		for( hd_size i=0; i<total_cand_count; ++i ) {
			hd_size bin = (hd_size)((h_all_dms[i] - dm_min) / dm_hist_binsep);
			++dm_hist[bin];
		}
		std::ofstream dm_hist_file("dm_hist.dat");
		for( hd_size b=0; b<dm_hist_nbins; ++b ) {
			dm_hist_file << (b+1)*dm_hist_binsep << "\t" << dm_hist[b] << "\n";
		}
		dm_hist_file.close();
		// -------------------------------------------
		*/


		/*
		for bi in beams:
	      for ci in bi.candidates:
	          for bj in beams:
	              if( bi == bj ) continue
	              for cj in bj.candidates:
	                  if( are_coincident(ci, cj) ):
	                      ci.beam_count += 1;
	                      ci.beam_mask |= 1 << bj;
	                      continue;
		*/


	/*
	  label_candidate_clusters(all events) --> labels
	  sort(values=zip(all events), keys=labels)
	  reduce_by_key(values=zip(all events), keys=labels) --> zip(event groups)
	  
	  Label event groups as NOISE if nmembers < 3
	    TODO: Check if this includes members across time samples
	            My understanding is that it doesn't
	  Label event groups as DM0 RFI if DMpeak < 1.5 pc cm^-3
	  
	  Somehow determine beam coincidence for each event group
	  Label event groups as COINC RFI if beam coincidence >= rfi_min_beams
	  
	  sorted = sort(all_candidates by time);
	  adjacent_difference(
	  
	  Send candidates to root node
	  On root node
	  ------------
	  for bi in beams:
	      for ci in bi.candidates:
	          for bj in beams:
	              if( bi == bj ) continue
	              for cj in bj.candidates:
	                  if( are_coincident(ci, cj) ):
	                      ci.beam_count += 1;
	                      ci.beam_mask |= 1 << bj;
	                      continue;
	  for ci in candidates:
	      if( ci.beam_coincidence > pl->params.rfi_min_beams ) {
	          ci.flags |= FLAG_MULTIBEAM_RFI;
	      }
	  
	  Gather together all events whose group is labelled as NOISE
	  Gather together all events whose group is labelled as DM0 RFI
	  Gather together all events whose group is labelled as COINC RFI
	  
	  Plot NOISE events, DM0 RFI events and COINC RFI events with different styles
	  
	  //Sort by SNRpeak (descending)
	  
	 */
	
// Repeat from (1) for each candidate data array


		// TODO: Get external commands and act accordingly
		//read_source_data(filterbank, &nsamps, &nbits);
		//DummyDataSource data_source("pulsar_dm0020.fil");
		//DummyDataSource data_source("transient_15s_8b.fil");
		
//int proc_idx;
		//MPI_Comm_rank(MPI_COMM_WORLD, &proc_idx);

//size_t first_beam_idx = 10 - 1;
		//size_t first_beam_idx = 11 - 1;

		//cout << "Beginning from beam " << params.first_beam + 1 << endl;

/*
// TODO: Remove this when satisfied with the OOP version below
// Note: This operates quite similarly to remove_baseline
// Note: This can operate "in-place"
hd_error matched_filter(const hd_float* d_in, hd_size count, hd_size filter_width,
                        hd_float* d_out)
{
	hd_size w = filter_width;
	
	thrust::device_ptr<const hd_float> d_in_begin(d_in);
	thrust::device_ptr<const hd_float> d_in_end(d_in + count);
	thrust::device_ptr<hd_float>       d_out_begin(d_out);
	thrust::device_ptr<hd_float>       d_out_end(d_out + count);
	
	// Check if there's actually work to do
	if( 1 == filter_width ) {
		if( d_in != d_out ) {
			thrust::copy(d_in_begin, d_in_end, d_out_begin);
		}
		return HD_NO_ERROR;
	}
	
	// TODO: Check if this repeated memory allocation is slow
	//         Could try to allocate it once only
	// Note: The ends are padded with zeros
	thrust::device_vector<hd_float> d_scanned(count + w, hd_float(0));
	
	// We copy the input into the middle of the padded d_scanned array
	thrust::copy(d_in_begin, d_in_end, d_scanned.begin()+w/2);
	thrust::exclusive_scan(d_scanned.begin(), d_scanned.end(),
	                       d_scanned.begin());
	
	thrust::transform(d_scanned.begin()+w, d_scanned.end(),
	                  d_scanned.begin(),
	                  d_out_begin,
	                  averager<hd_float>(w));
	
	return HD_NO_ERROR;
}
*/

