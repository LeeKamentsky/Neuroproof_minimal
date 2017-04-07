#include "vigraRFclassifier.h"
// #include "assert.h"
#include <time.h>
#include <vector>

namespace vigra {


// NOTE(TFK): Highly conservative threshold check. Improves agglomeration time by a factor of 2.
     /// should never change the result.
class StopConservativeThreshold : public StopBase
{
public:
    double max_tree_p;
    int max_tree_;
    typedef StopBase SB;
    int min_split_node_size_;
    double thresh;
    int tree_count;

    template<class Opt>
    StopConservativeThreshold(Opt opt, double _thresh) {
      thresh = _thresh;
      min_split_node_size_ = opt.min_split_node_size_; 
      tree_count = opt.tree_count_; 
    }

    template<class Region>
    bool operator()(Region& region)
    {
        return region.size() < min_split_node_size_;
    }

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter weightIter, int k, MultiArrayView<2, T, C> const &  prob , double totalCt) {
      if (prob(0,1) > thresh) return true;
      return false;
    }
};
}

VigraRFclassifier::VigraRFclassifier(const char* rf_filename){

    _rf=NULL;	
    load_classifier(rf_filename);
}

void  VigraRFclassifier::load_classifier(const char* rf_filename){
      HDF5File rf_file(rf_filename,HDF5File::Open); 	
      if(!_rf)
	      _rf = new RandomForest<>;
      rf_import_HDF5(*_rf, rf_file,"rf");
      printf("RF loaded with %d trees, for %d class prob with %d dimensions\n",_rf->tree_count(), _rf->class_count(), _rf->column_count());
      _nfeatures = _rf->column_count();
      _nclass = _rf->class_count();
      

      /* read list of useless features*/ 
      string filename = rf_filename;
      unsigned found = filename.find_last_of(".");
      string nameonly = filename.substr(0,found);	
      nameonly += "_ignore.txt";
      FILE* fp = fopen(nameonly.c_str(),"rt");
      if (!fp){
	  printf("no features to ignore\n");
	  return;
      }
      else{
	unsigned int a;
	while(!feof(fp)){
	    fscanf(fp,"%u ", &a);
	    ignore_featlist.push_back(a);
	    
	}
	fclose(fp);
      }
  
	
}

std::vector<double> VigraRFclassifier::predict_batch(std::vector< std::vector<double> >& batch_features){
	//if(!_rf){
	//    return(EdgeClassifier::predict(features));
	//}
// 	if(_nfeatures != features.size()){
// 	    printf("number of features differet in classifer\n ");
// 	    exit(0);
// 	}

        MultiArray<2, float> vfeatures(Shape(batch_features.size(),_nfeatures));
     	MultiArray<2, float> prob(Shape(batch_features.size(), _nclass));

        for (int j = 0; j < batch_features.size(); j++) {
	  for(int i = 0; i < _nfeatures; i++) {
	        vfeatures(j,i)= (float)batch_features[j][i];
          }
        }
        vigra::StopConservativeThreshold stop = vigra::StopConservativeThreshold(_rf->options_,0.2*255);
        _rf->predictProbabilities(vfeatures, prob, stop);

	/*debug*/
	std::vector<double> vp(_nclass);
	for(int i=0;i<_nclass;i++)
	    vp[i] = prob(0,i);	
	/**/
        std::vector<double> results;
        for (int i = 0; i < batch_features.size(); i++) {
          results.push_back(prob(i,1));
        }
        return results;
	//return (double) prob(0,1);
}	

double VigraRFclassifier::predict(std::vector<double>& features){
	if(!_rf){
	    return(EdgeClassifier::predict(features));
	}
// 	if(_nfeatures != features.size()){
// 	    printf("number of features differet in classifer\n ");
// 	    exit(0);
// 	}
        MultiArray<2, float> vfeatures(Shape(1,_nfeatures));
     	MultiArray<2, float> prob(Shape(1, _nclass));
	for(int i=0;i<_nfeatures;i++)
	      vfeatures(0,i)= (float)features[i];

        vigra::StopConservativeThreshold stop = vigra::StopConservativeThreshold(_rf->options_, 0.2*255);
        _rf->predictProbabilities(vfeatures, prob, stop);
        //_rf->predictProbabilities(vfeatures, prob);    

	/*debug*/
	std::vector<double> vp(_nclass);
	for(int i=0;i<_nclass;i++)
	    vp[i] = prob(0,i);	
	/**/

	return (double) prob(0,1);
			    	
}	


void VigraRFclassifier::learn(std::vector< std::vector<double> >& pfeatures, std::vector<int>& plabels){

     if (_rf)
	delete _rf;	

     unsigned int rows = pfeatures.size();
     unsigned int cols = pfeatures[0].size();	 	
     
     printf("Number of samples and dimensions: %u, %u\n",rows, cols);
     if ((rows<1)||(cols<1)){
	return;
     }

     clock_t start = clock();

     MultiArray<2, float> features(Shape(rows,cols));
     MultiArray<2, int> labels(Shape(rows,1));

     unsigned int numzeros=0; 	
     for(unsigned int i=0; i < rows; i++){
	 labels(i,0) = plabels[i];
	 numzeros += (labels(i,0)==-1?1:0);
	 for(unsigned int j=0; j < cols ; j++){
	     features(i,j) = (float)pfeatures[i][j];	
	 }
     }
     printf("Number of merge: %u\n",numzeros);



     int tre_count = 255;	
     printf("Number of trees:  %u\n",tre_count);
     RandomForestOptions rfoptions = RandomForestOptions().tree_count(tre_count).use_stratification(RF_EQUAL);	//RF_EQUAL, RF_PROPORTIONAL
     _rf= new RandomForest<>(rfoptions);


     // construct visitor to calculate out-of-bag error
     visitors::OOB_Error oob_v;
     visitors::VariableImportanceVisitor varimp_v;

     _rf->learn(features, labels);
     _nfeatures = _rf->feature_count();
     _nclass = _rf->class_count();

     printf("Time required to learn RF: %.2f with oob :%f\n", ((double)clock() - start) / CLOCKS_PER_SEC, oob_v.oob_breiman);
}

void VigraRFclassifier::save_classifier(const char* rf_filename){

  
    if (ignore_featlist.size()>0){
	string filename = rf_filename;
	unsigned found = filename.find_last_of(".");
	string nameonly = filename.substr(0,found);	
	nameonly += "_ignore.txt";
	FILE* fp = fopen(nameonly.c_str(),"wt");
	for(size_t i=0; i<ignore_featlist.size(); i++){
	    fprintf(fp,"%u ", ignore_featlist[i]);
	}
	fclose(fp);
    }
  
  
     HDF5File rf_file(rf_filename,HDF5File::New); 	
     rf_export_HDF5(*_rf, rf_file,"rf");
}

