#include "DataStructures/Stack.h"
// #include "Priority/GPR.h"
// #include "Priority/LocalEdgePriority.h"
// #include "Utilities/ScopeTime.h"
// #include "ImportsExports/ImportExportRagPriority.h"
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <memory>
// #include <json/json.h>
// #include <json/value.h>

#include "Utilities/h5read.h"
#include "Utilities/h5write.h"

#include <ctime>
#include <cmath>
#include <cstring>

//#include <google/profiler.h>

using std::cerr; using std::cout; using std::endl;
using std::ifstream;
using std::string;
using std::stringstream;
using namespace NeuroProof;
using std::vector;


extern float array_items[256];

template <class T>
void padZero(T* data, const size_t* dims, int pad_length, T** ppadded_data){
     
     // implemented only for 3D arrays

     unsigned long int newsize = (dims[0]+2*pad_length)*(dims[1]+2*pad_length)*(dims[2]+2*pad_length);	
     *ppadded_data = new T[newsize];  	
     T* padded_data = *ppadded_data;
	
     memset((void*) padded_data, 0, sizeof(T)*newsize);

     unsigned int width, plane_size, width0, plane_size0, i0,j0,k0, i,j,k;
     
     for (i=pad_length, i0=0; i<= dims[0] ; i++, i0++)
	 for (j=pad_length, j0=0; j<= dims[1]; j++, j0++)	
	     for(k=pad_length, k0=0; k<= dims[2]; k++, k0++){
		 plane_size = (dims[1]+2*pad_length)*(dims[2]+2*pad_length);	
		 width = dims[2]+2*pad_length;	

		 plane_size0 = dims[1]*dims[2];	
		 width0 = dims[2];	

		 padded_data[i*plane_size+ j*width + k] = data[i0*plane_size0 + j0*width0 + k0];	
	     }
			

}

bool endswith(string filename, string extn){
    unsigned found = filename.find_last_of(".");
    string fextn = filename.substr(found);	
    if (fextn.compare(extn) == 0 )
	return true;
    else return false;	  
}


EdgeClassifier* load_classifier(std::string classifier_filename) {
    EdgeClassifier* eclfr;
    if (endswith(classifier_filename, ".h5")){
	string nameonly = classifier_filename.substr(0, classifier_filename.find_last_of("."));	
	    eclfr = new VigraRFclassifier(classifier_filename.c_str());	
    }	      
    else if (endswith(classifier_filename, ".xml")){
	string nameonly = classifier_filename.substr(0, classifier_filename.find_last_of("."));	
	    eclfr = new OpencvRFclassifier(classifier_filename.c_str());	
    }
    return eclfr;
}

void load_watershed(std::string watershed_filename, std::string watershed_dataset_name, Label** zp_watershed_data) {
    H5Read watershed(watershed_filename.c_str(),watershed_dataset_name.c_str());	
    Label* watershed_data=NULL;	
    watershed.readData(&watershed_data);	
    int depth =	 watershed.dim()[0];
    int height = watershed.dim()[1];
    int width =	 watershed.dim()[2];
	
    int pad_len=1;
    //Label *zp_watershed_data=NULL;
    padZero(watershed_data, watershed.dim(),pad_len,zp_watershed_data);

    // NOTE(TFK): Original unpadded watershed data is never used after this point. Can delete.
    delete[] watershed_data;

}



int main(int argc, char** argv) 
{
    int          i, j, k;


    cout<< "Reading data ..." <<endl;
    printf("%f\n", 1.0/256);


    if (argc<15){
	printf("format: NeuroProof_stack -watershed watershed_h5_file  dataset \
					 -prediction prediction_h5_file dataset \
					 -classifier classifier_file \
					 -output output_h5_file dataset \
					 -groundtruth groundtruth_h5_file dataset \
					 -threshold agglomeration_threshold \
					 -algorithm algorithm_type\n");
	return 0;
    }	

    //ProfilerStart("profile.data");
    int argc_itr=1;	
    string watershed_filename="";
    string watershed_dataset_name="";		
    string prediction_filename="";
    string prediction_dataset_name="";		
    string output_filename;
    string output_dataset_name;		
    string groundtruth_filename="";
    string groundtruth_dataset_name="";		
    string classifier_filename;

    string output_filename_nomito;

    	
    double threshold = 0.2;	
    int agglo_type = 1;		
    bool merge_mito = true;
    bool merge_mito_by_chull = false;
    bool read_off_rwts = false;
    double mito_thd=0.3;
    size_t min_region_sz=100;
    while(argc_itr<argc){
	if (!(strcmp(argv[argc_itr],"-watershed"))){
	    watershed_filename = argv[++argc_itr];
	    watershed_dataset_name = argv[++argc_itr];
	}

	if (!(strcmp(argv[argc_itr],"-classifier"))){
	    classifier_filename = argv[++argc_itr];
	}
	if (!(strcmp(argv[argc_itr],"-prediction"))){
	    prediction_filename = argv[++argc_itr];
	    prediction_dataset_name = argv[++argc_itr];
	}
	if (!(strcmp(argv[argc_itr],"-output"))){
	    output_filename = argv[++argc_itr];
	    output_dataset_name = argv[++argc_itr];
	}
	if (!(strcmp(argv[argc_itr],"-groundtruth"))){
	    groundtruth_filename = argv[++argc_itr];
	    groundtruth_dataset_name = argv[++argc_itr];
	}

	if (!(strcmp(argv[argc_itr],"-threshold"))){
	    threshold = atof(argv[++argc_itr]);
	}
	if (!(strcmp(argv[argc_itr],"-algorithm"))){
	    agglo_type = atoi(argv[++argc_itr]);
	}
	if (!(strcmp(argv[argc_itr],"-min_region_sz"))){
	    min_region_sz = atoi(argv[++argc_itr]);
	}
	if (!(strcmp(argv[argc_itr],"-nomito"))){
	    merge_mito = false;
	}
	if (!(strcmp(argv[argc_itr],"-mito_thd"))){
	    mito_thd = atof(argv[++argc_itr]);
	}
	if (!(strcmp(argv[argc_itr],"-read_off"))){
	    if (agglo_type==2)
		read_off_rwts = true;
	}
        ++argc_itr;
    } 	
    	
    output_filename_nomito = output_filename;
    size_t ofn = output_filename_nomito.find_last_of(".");
    output_filename_nomito.erase(ofn-1,1);			

    time_t start, end, start_agg, start_rag;
    time(&start);	

    EdgeClassifier* eclfr = load_classifier(classifier_filename);
//    if (endswith(classifier_filename, ".h5")){
//	string nameonly = classifier_filename.substr(0, classifier_filename.find_last_of("."));	
//// 	if (nameonly.find("parallel") != std::string::npos)
//// 	    eclfr = new VigraRFclassifierP(classifier_filename.c_str());
//// 	else
//	    eclfr = new VigraRFclassifier(classifier_filename.c_str());	
//    }	      
////     	eclfr = new VigraRFclassifier(classifier_filename.c_str());	
//    else if (endswith(classifier_filename, ".xml")){
//	string nameonly = classifier_filename.substr(0, classifier_filename.find_last_of("."));	
//// 	if (nameonly.find("parallel") != std::string::npos)
//// 	    eclfr = new OpencvRFclassifierP(classifier_filename.c_str());	
//// 	else
//	    eclfr = new OpencvRFclassifier(classifier_filename.c_str());	
//    }
//	//eclfr = new CompositeRFclassifier(classifier_filename.c_str());	
//// 	eclfr = new OpencvRFclassifier(classifier_filename.c_str());	

    

    H5Read prediction(prediction_filename.c_str(),prediction_dataset_name.c_str(), true);	
    float* prediction_data=NULL;
    prediction.readData(&prediction_data);	
    int depth =	 prediction.dim()[0];
    int height = prediction.dim()[1];
    int width =	 prediction.dim()[2];
    int pad_len=1;
	
    //float* prediction_ch0 = new float[depth*height*width];


    // NOTE(TFK): prediction_single_ch is used as a buffer. The last channel (ch0) is
    //    process last, and so prediction_ch0 points to it. This is correct because
    //    after processing all channels the buffer will contain the contents of channel 0.
    unsigned char* prediction_single_ch = new unsigned char[depth*height*width];
    unsigned char* prediction_ch0 = prediction_single_ch;

    size_t cube_size, plane_size, element_size=prediction.dim()[prediction.total_dim()-1];
    size_t position, count;		  


    // NOTE(TFK): In order to reuse the buffer space (prediction_single_ch) for prediction_ch0
    //    we need to process the channels in reverse order. To avoid changing the
    //    order of prediction channels in the Stack we push the channels to a vector
    //    and add them to the Stack in order of channel id.  	
    //for (int ch=0; ch < element_size; ch++){
    std::vector<unsigned char*> prediction_channel_list;
    //std::set<float> dedupe;
    for (int ch = element_size; --ch >=0;) {
	count = 0;
	for(int i=0; i<depth; i++){
	    int cube_size = height*width*element_size;	
	    for(int j=0; j<height; j++){
		int plane_size = width*element_size;
		for(int k=0; k<width; k++){
		    int position = i*cube_size + j*plane_size + k*element_size + ch;
                    if (ch==element_size-1) {
		      //prediction_single_ch[count] = (unsigned char) prediction_data[position]; 
                      array_items[(unsigned char) lround(255-255*prediction_data[position-1])] = prediction_data[position];
                      //if (dedupe.find(prediction_data[position]) == dedupe.end()) {
                      //  dedupe.insert(prediction_data[position]);
                      //  printf("entry %d value %f\n", lround(255-255*prediction_data[position-1]), prediction_data[position]);
                      //}
                    } else {
                      //int real_post = i*height*width + j*width + k;
		      prediction_single_ch[count] = (unsigned char) lround(255*prediction_data[position]);
                    }
		    count++;					
                    //__sync_fetch_and_add(&count, 1);
		}		
	    }	
	}

	unsigned char* zp_prediction_single_ch = NULL;
        if (ch==0) {	
	  padZero(prediction_single_ch,prediction.dim(),pad_len,&zp_prediction_single_ch);
        } 
          prediction_channel_list.push_back(zp_prediction_single_ch);
          
	//if (ch == 0)
	//    memcpy(prediction_ch0, prediction_single_ch, depth*height*width*sizeof(float));	
	//stackp->add_prediction_channel(zp_prediction_single_ch);
    }

    //printf("size fo the set is %d\n", dedupe.size());

    // NOTE(TFK): Prediction data not used after this point. This is not ideal place to
    //     delete this data, because I don't believe it decreases the code's memory
    //     high-water mark. But might as well delete it.
    delete[] prediction_data;

    Label* zp_watershed_data = NULL;
    load_watershed(watershed_filename, watershed_dataset_name, &zp_watershed_data);
    StackPredict* stackp = new StackPredict(zp_watershed_data, depth+2*pad_len, height+2*pad_len, width+2*pad_len, pad_len);
    stackp->set_feature_mgr(new FeatureMgr());
    stackp->set_merge_mito(merge_mito, mito_thd);

    for (int i = prediction_channel_list.size(); --i>=0;) {
      stackp->add_prediction_channel(prediction_channel_list[i]);
    }

    stackp->set_basic_features();
    stackp->get_feature_mgr()->set_classifier(eclfr);   	 

    Label* groundtruth_data=NULL;
    if (groundtruth_filename.size()>0){  	
	H5Read groundtruth(groundtruth_filename.c_str(),groundtruth_dataset_name.c_str());	
        groundtruth.readData(&groundtruth_data);	
        stackp->set_groundtruth(groundtruth_data);
    }	
	

    cout<<"Building RAG ..."; 	
    time(&start_rag);	
    stackp->build_rag();
    cout<<"done with "<< stackp->get_num_bodies()<< " regions\n";	
    cout<<"Inclusion removal ..."; 
    time(&end);	
    printf("Build Rag Time: %.4f\n", (difftime(end,start_rag))*1.0);
    stackp->remove_inclusions();
    cout<<"done with "<< stackp->get_num_bodies()<< " regions\n";	

    time(&start_rag);	
    stackp->compute_vi();  	
    stackp->compute_groundtruth_assignment();
    time(&end);	
    printf("Compute vi and groundtruth assignment time: %.4f\n", (difftime(end,start_rag))*1.0);
    time(&start_agg);	
    if (agglo_type==0){	
	cout<<"Agglomerating (flat) upto threshold "<< threshold<< " ...\n"; 
    	stackp->agglomerate_rag_flat(threshold);	
    }
    else if (agglo_type==1){
	cout<<"Agglomerating (agglo) upto threshold "<< threshold<< " ...\n"; 
     	stackp->agglomerate_rag(threshold, false);	
    }
    else if (agglo_type == 2){		
	cout<<"Agglomerating (mrf) upto threshold "<< threshold<< " ...\n"; 
    	stackp->agglomerate_rag_mrf(threshold, read_off_rwts, output_filename, classifier_filename);	
    }		
    else if (agglo_type == 3){		
	cout<<"Agglomerating (queue) upto threshold "<< threshold<< " ...\n"; 
    	stackp->agglomerate_rag_queue(threshold, false);	
    }		
    else if (agglo_type == 4){		
	cout<<"Agglomerating (flat) upto threshold "<< threshold<< " ...\n"; 
    	stackp->agglomerate_rag_flat(threshold, false, output_filename, classifier_filename);	
    }		
    cout << "Done with "<< stackp->get_num_bodies()<< " regions\n";
    time(&end);	
    printf("Agglomeration Time: %.4f\n", (difftime(end,start_agg))*1.0);

    time(&start_agg);	
    cout<<"Inclusion removal ..."; 
    stackp->remove_inclusions();
    cout<<"done with "<< stackp->get_num_bodies()<< " regions\n";	
    time(&end);	
    printf("Remove inclusion time: %.4f\n", (difftime(end,start_agg))*1.0);



    stackp->compute_vi();  	
 	
    time(&start_agg);	
    Label * temp_label_volume1D = stackp->get_label_volume();       	    
    if (!merge_mito/* && min_region_sz > 0*/) {
      printf("Absort small regions2 in !merge_mito\n");
      //stackp->absorb_small_regions2(prediction_ch0, temp_label_volume1D, min_region_sz);
      printf("After Absort small regions2 in !merge_mito\n");
    }

    hsize_t dims_out[3];
    dims_out[0]=depth; dims_out[1]= height; dims_out[2]= width;   
    H5Write(output_filename_nomito.c_str(),output_dataset_name.c_str(),3,dims_out, temp_label_volume1D);
    printf("Output-nomito written to %s, dataset %s\n",output_filename_nomito.c_str(),output_dataset_name.c_str());	
    delete[] temp_label_volume1D;	
    time(&end);	
    printf("Writing output time %.4f\n", (difftime(end,start_agg))*1.0);

    if (merge_mito){
	cout<<"Merge Mitochondria (border-len) ..."; 
	stackp->merge_mitochondria_a();
    	cout<<"done with "<< stackp->get_num_bodies()<< " regions\n";	

        cout<<"Inclusion removal ..."; 
        stackp->remove_inclusions();
        cout<<"done with "<< stackp->get_num_bodies()<< " regions\n";	

        stackp->compute_vi();  	
 	
        temp_label_volume1D = stackp->get_label_volume();       	    
        printf("Absorb small regions2 in merge_mito\n");
        // NOTE(TFK): stackp->absorb_small_regions2(prediction_ch0, temp_label_volume1D, min_region_sz);
        printf("After Absorb small regions2 in merge_mito\n");

        dims_out[0]=depth; dims_out[1]= height; dims_out[2]= width;   
        H5Write(output_filename.c_str(),output_dataset_name.c_str(),3,dims_out, temp_label_volume1D);
        printf("Output written to %s, dataset %s\n",output_filename.c_str(),output_dataset_name.c_str());	
        delete[] temp_label_volume1D;	
    } 	
    //stackp->determine_edge_locations();
    //stackp->write_graph(output_filename);

    time(&end);	
    printf("Time elapsed: %.2f\n", (difftime(end,start))*1.0/60);



//    if (watershed_data)  	
//	delete[] watershed_data;
    delete[] zp_watershed_data;
	
    //if (prediction_data)  	
    //  delete[] prediction_data;
    delete[] prediction_single_ch;

    if (groundtruth_data)  	
	delete[] groundtruth_data;

     	
    //ProfilerStop();

    return 0;
}
