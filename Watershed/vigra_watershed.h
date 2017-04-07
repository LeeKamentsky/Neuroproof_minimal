#ifndef _vigra_watershed
#define _vigra_watershed


#include <vigra/seededregiongrowing3d.hxx>
#include <vigra/multi_array.hxx>

using namespace vigra;
 

typedef vigra::MultiArray<3,unsigned int> LblVolume;
typedef vigra::MultiArray<3,double> DVolume;


class VigraWatershed{
    size_t _depth; 	
    size_t _width; 	
    size_t _height; 	

public:
    VigraWatershed(size_t depth, size_t height, size_t width): _depth(depth), _height(height), _width(width) {};  	

    #define VIGRA_WATER false
    void run_watershed(float* data_vol, unsigned int* lbl_vol){
        #if VIGRA_WATER
    	DVolume src(DVolume::difference_type(_depth, _height, _width));
	LblVolume dest(LblVolume::difference_type(_depth, _height, _width));
        #endif

    	unsigned long volsz = _depth*_height*_width;


        unsigned int width = _width;

    	unsigned long nnz=0; 	
    	
/*
	// *C* for debug only
	for(unsigned long i=0; i< volsz; i++){
	    unsigned int lbl = lbl_vol[i];
	    if (lbl==0){
	    	nnz++;	
	    }
    	} 
    	printf("total #zeros 1D: %lu\n",nnz);
*/

	unsigned long plane_size = _height*_width;
	nnz = 0;

        std::priority_queue<std::pair<float,unsigned long> > seed_queue;
        std::set<unsigned long> seed_queue_set;
        unsigned long max_label = 0;
	for(unsigned int d=0; d < _depth; d++){
	    unsigned long z1d = d*plane_size;
	    for(unsigned int i=0; i < _height; i++){
		unsigned long x1d = i*_width;
		for(unsigned int j=0; j < _width; j++){
		    unsigned long l1d = z1d+ x1d + j;
                    #if VIGRA_WATER
		    src(d,i,j) = data_vol[l1d];
		    dest(d,i,j) = lbl_vol[l1d];		
                    #endif	
                    if (lbl_vol[l1d] > max_label) max_label = lbl_vol[l1d];

		    if (lbl_vol[l1d] == 0){
		    //if (dest(d,i,j) == 0){
			int pp=1;
			nnz++;

                        // Find seeds adjacent to a region.
                        bool adj = false;

                        // InFront
                        //if (l1d >= plane_size) {
                        //  if (lbl_vol[l1d-plane_size] != 0) {
                        //    adj = true;
                        //  } 
                        //}

                        //if (l1d + plane_size < volsz) {
                        //  if (lbl_vol[l1d+plane_size] != 0) {
                        //    adj = true;
                        //  } 
                        //}

                        //if (l1d >= width) {
                        //  if (lbl_vol[l1d-width] != 0) {
                        //    adj = true;
                        //  } 
                        //}
                        //if (l1d + width < volsz) {
                        //  if (lbl_vol[l1d+width] != 0) {
                        //    adj = true;
                        //  } 
                        //}
//                        if (l1d >= 1) {
//                          if (lbl_vol[l1d-1] != 0) {
//                            adj = true;
//                          } 
//                        }


                        if (l1d + 1 < volsz) {
                          if (lbl_vol[l1d+1] != 0) {
                            adj = true;
                          } 
                        }
                        if (adj) {
                          seed_queue.push(std::make_pair(-1.0*data_vol[l1d], l1d));
                          seed_queue_set.insert(l1d);
                          //seed_list.push_back(l1d);	
                        }
		    }
		}
	    }	 	
	}


        while (!(VIGRA_WATER)&&!seed_queue.empty()) {
          // index into lbl_vol and data_vol.
          unsigned long l1d = seed_queue.top().second;
          seed_queue.pop();



          // Behind
          //if (l1d + plane_size < volsz) {
          //  if (lbl_vol[l1d+plane_size] != 0) {
          //    lbl_vol[l1d] = lbl_vol[l1d+plane_size];
          //    goto update; 
          //  } 
          //}

          // InFront
          //if (l1d >= plane_size) {
          //  if (lbl_vol[l1d-plane_size] != 0) {
          //    lbl_vol[l1d] = lbl_vol[l1d-plane_size];
          //    goto update;
          //  } 
          //}
          // North
          //if (l1d >= width) {
          //  if (lbl_vol[l1d-width] != 0) {
          //    lbl_vol[l1d] = lbl_vol[l1d-width];
          //    goto update; 
          //  } 
          //}          // South
          //if (l1d + width < volsz) {
          //  if (lbl_vol[l1d+width] != 0) {
          //    lbl_vol[l1d] = lbl_vol[l1d+width];
          //    goto update; 
          //  } 
          //}

          // West
          //if (l1d >= 1) {
          //  if (lbl_vol[l1d-1] != 0) {
          //    lbl_vol[l1d] = lbl_vol[l1d-1];
          //    goto update; 
          //  } 
          //}








          // East
          if (l1d + 1 < volsz) {
            if (lbl_vol[l1d+1] != 0) {
              lbl_vol[l1d] = lbl_vol[l1d+1];
              goto update; 
            } 
          }





        update:

//                        if (l1d >= width) {
//                          if (lbl_vol[l1d-width] == 0) {
//                            if (seed_queue_set.find(l1d-width) == seed_queue_set.end()) {
//                              seed_queue_set.insert(l1d-width);
//                              seed_queue.push(std::make_pair(-1.0*data_vol[l1d-width], l1d-width));
//                            } 
//                          } 
//                        }
//                        if (l1d + width < volsz) {
//                          if (lbl_vol[l1d+width] == 0) {
//                            if (seed_queue_set.find(l1d+width) == seed_queue_set.end()) {
//                              seed_queue_set.insert(l1d+width);
//                              seed_queue.push(std::make_pair(-1.0*data_vol[l1d+width], l1d+width));
//                            } 
//                          } 
//                        }
//
//                        if (l1d >= 1) {
//                          if (lbl_vol[l1d-1] == 0) {
//                            if (seed_queue_set.find(l1d-1) == seed_queue_set.end()) {
//                              seed_queue_set.insert(l1d-1);
//                              seed_queue.push(std::make_pair(-1.0*data_vol[l1d-1], l1d-1));
//                            } 
//                          } 
//                        }

                        if (l1d + 1 < volsz) {
                          if (lbl_vol[l1d+1] == 0) {
                            if (seed_queue_set.find(l1d+1) == seed_queue_set.end()) {
                              seed_queue_set.insert(l1d+1);
                              seed_queue.push(std::make_pair(-1.0*data_vol[l1d+1], l1d+1));
                            } 
                          } 
                        }
			//if (l1d + plane_size < volsz) {
                        //  if (lbl_vol[l1d+plane_size] == 0) {
                        //    if (seed_queue_set.find(l1d+plane_size) == seed_queue_set.end()) {
                        //      seed_queue_set.insert(l1d+plane_size);
                        //      seed_queue.push(std::make_pair(-1.0*data_vol[l1d+plane_size], l1d+plane_size));
                        //    } 
                        //  } 
                        //}

                        // InFront
                        //if (l1d >= plane_size) {
                        //  if (lbl_vol[l1d-plane_size] == 0) {
                        //    if (seed_queue_set.find(l1d-plane_size) == seed_queue_set.end()) {
                        //      seed_queue_set.insert(l1d-plane_size);
                        //      seed_queue.push(std::make_pair(-1.0*data_vol[l1d-plane_size], l1d-plane_size));
                        //    } 
                        //  } 
			//}
/*

			if (l1d + plane_size < volsz) {
                          if (lbl_vol[l1d+plane_size] == 0) {
                            if (seed_queue_set.find(l1d+plane_size) == seed_queue_set.end()) {
                              seed_queue_set.insert(l1d+plane_size);
                              seed_queue.push(std::make_pair(-1.0*data_vol[l1d+plane_size], l1d+plane_size));
                            } 
                          } 
                        }
*/
        }
        

        
#if VIGRA_WATER
    	vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<double> >
                                              stats;

    	
    vigra::seededRegionGrowing3D(srcMultiArrayRange(src), destMultiArray(dest),
                               destMultiArray(dest), stats);
#endif

	unsigned long nnz2=0;
	for(size_t d=0; d < _depth; d++){
	    unsigned long z1d = d*plane_size;
	    for(size_t i=0; i < _height; i++){
		unsigned long x1d = i*_width;
		for(size_t j=0; j < _width; j++){
		    unsigned long l1d = z1d+ x1d + j;
                    #if VIGRA_WATER
		    lbl_vol[l1d] = dest(d,i,j);		
                    #endif	
		    if (lbl_vol[l1d] == 0){
			int pp=1;
			nnz2++;
		    }
		}	
	    }	 	
	}
    	printf("total #zeros in vol before and after wshed: %lu,  %lu\n",nnz, nnz2);
    }	
};

#endif
