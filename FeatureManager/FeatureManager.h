#ifndef FEATUREMANAGER_H
#define FEATUREMANAGER_H


#include "../DataStructures/RagEdge.h"
#include "../DataStructures/Glb.h"
#include "Features.h"
#include <tr1/unordered_map>

#define SIMULATE_SECOND_CHANNEL true

// #ifdef SETPYTHON
// #include <boost/python.hpp>
// #endif


#include "../Classifier/vigraRFclassifier.h"
// #include "../Classifier/vigraRFclassifierP.h"
#include "../Classifier/opencvRFclassifier.h"
// #include "../Classifier/opencvRFclassifierP.h"
// #include "../Classifier/compositeRFclassifier.h"
// #include "../Classifier/opencvSVMclassifier.h"

// node features and edge features -- different amounts ??!!


namespace NeuroProof {

typedef std::tr1::unordered_map<RagEdge<Label>*, std::vector<void *>, RagEdgePtrHash<Label>, RagEdgePtrEq<Label> > EdgeCaches; 
typedef std::tr1::unordered_map<RagNode<Label>*, std::vector<void *>, RagNodePtrHash<Label>, RagNodePtrEq<Label> > NodeCaches; 


class FeatureMgr {
  public:
    FeatureMgr() : num_channels(0), specified_features(false), has_pyfunc(false), num_features(0), eclfr(NULL) {}
    FeatureMgr(int num_channels_) : num_channels(num_channels_), specified_features(false), channels_features(num_channels_), channels_features_modes(num_channels_), has_pyfunc(false), num_features(0), eclfr(NULL) {}
    
    void add_channel();
    unsigned int get_num_features()
    {
        return num_features;
    }
#ifdef SETPYTHON
    void set_python_rf_function(boost::python::object pyfunc_);
#endif

    unsigned int get_num_channels()
    {
        return num_channels;
    }

    void set_simulate(bool do_simulate) {
	simulate = do_simulate;
    }
    bool simulate_channel_1() {
	return simulate;
    }
    
    void add_val(unsigned char _val, RagNode<Label>* node)
    {
        float val = _val*1.0/255.0;
        unsigned int starting_pos = 0;
        if (node_caches.find(node) != node_caches.end()) {
            add_val(val, 0, starting_pos, node_caches[node]);
        } else {
            std::vector<void*>& feature_caches = create_cache(node);
            add_val(val, 0, starting_pos, feature_caches);
        } 
    } 
    
    void add_val(unsigned char _val, RagEdge<Label>* edge)
    {
        float val = _val*1.0/255.0;
        unsigned int starting_pos = 0;
        if (edge_caches.find(edge) != edge_caches.end()) {
            add_val(val, 0, starting_pos, edge_caches[edge]);
        } else {
            std::vector<void*>& feature_caches = create_cache(edge);
            add_val(val, 0, starting_pos, feature_caches);
        } 
    } 
 
    void add_val_batch(std::vector<unsigned char>* val, unsigned int channel, unsigned int& starting_pos, std::vector<void *>& feature_caches)
    {
        std::vector<FeatureCompute*>& features = channels_features[channel];
        for (int i = 0; i < features.size(); ++i) {
	    if (simulate) {
		features[i]->add_point_batch(val, feature_caches[starting_pos], channel > 0); 
	    } else {
		features[i]->add_point_batch(val, feature_caches[starting_pos]);
	    }
            ++starting_pos;
        }
    }

    // for simulate second channel.
    void add_val_batch(std::vector<unsigned char>* values, RagNode<Label>* node)
    {
        node->incr_size(values->size());
        //assert(values->size() == num_channels);
        unsigned starting_pos = 0;
        if (node_caches.find(node) != node_caches.end()) {
            std::vector<void*>& feature_caches = node_caches[node];
            for (int i = 0; i < num_channels; ++i) { 
                add_val_batch(values, i, starting_pos, feature_caches);
            }              
        } else {
            std::vector<void*>& feature_caches = create_cache(node);
            for (int i = 0; i < num_channels; ++i) { 
                add_val_batch(values, i, starting_pos, feature_caches);
            }
        }        
    }

    void add_val_batch(std::vector<unsigned char>* values, RagEdge<Label>* edge)
    { 
        edge->incr_size(values->size());
        //assert(vals.size() == num_channels);
        unsigned int starting_pos = 0;
        if (edge_caches.find(edge) != edge_caches.end()) {
            std::vector<void*>& feature_caches = edge_caches[edge];
            for (int i = 0; i < num_channels; ++i) { 
                add_val_batch(values, i, starting_pos, feature_caches);
            }           
        } else {
            std::vector<void*>& feature_caches = create_cache(edge);
            for (int i = 0; i < num_channels; ++i) { 
                add_val_batch(values, i, starting_pos, feature_caches);
            }        
        } 
    }



    void add_val_batch(std::vector<std::vector<unsigned char> >* values, RagNode<Label>* node)
    {
        node->incr_size((*values)[0].size());
        assert(values->size() == num_channels);
        unsigned starting_pos = 0;
        if (node_caches.find(node) != node_caches.end()) {
            std::vector<void*>& feature_caches = node_caches[node];
            for (int i = 0; i < num_channels; ++i) { 
                add_val_batch(&((*values)[i]), i, starting_pos, feature_caches);
            }              
        } else {
            std::vector<void*>& feature_caches = create_cache(node);
            for (int i = 0; i < num_channels; ++i) { 
                add_val_batch(&((*values)[i]), i, starting_pos, feature_caches);
            }
        }        
    }

    void add_val_batch(std::vector<std::vector<unsigned char> >* values, RagEdge<Label>* edge)
    {
	edge->incr_size((*values)[0].size());
	assert(values->size() == num_channels);
	unsigned starting_pos = 0;
	if (edge_caches.find(edge) != edge_caches.end()) {
	    std::vector<void*>& feature_caches = edge_caches[edge];
	    for (int i = 0; i < num_channels; ++i) { 
		add_val_batch(&((*values)[i]), i, starting_pos, feature_caches);
	    }              
	} else {
	    std::vector<void*>& feature_caches = create_cache(edge);
	    for (int i = 0; i < num_channels; ++i) { 
		add_val_batch(&((*values)[i]), i, starting_pos, feature_caches);
	    }
	}        
    }

    void add_val(std::vector<unsigned char>& vals, RagNode<Label>* node)
    {
        node->incr_size();
        assert(vals.size() == num_channels);
        unsigned starting_pos = 0;
        if (node_caches.find(node) != node_caches.end()) {
            std::vector<void*>& feature_caches = node_caches[node];
            for (int i = 0; i < num_channels; ++i) { 
                add_val(vals[i], i, starting_pos, feature_caches);
            }              
        } else {
            std::vector<void*>& feature_caches = create_cache(node);
            for (int i = 0; i < num_channels; ++i) { 
                add_val(vals[i], i, starting_pos, feature_caches);
            }
        }        
    }

    void add_val(std::vector<unsigned char>& vals, RagEdge<Label>* edge)
    { 
        edge->incr_size();
        assert(vals.size() == num_channels);
        unsigned int starting_pos = 0;
        if (edge_caches.find(edge) != edge_caches.end()) {
            std::vector<void*>& feature_caches = edge_caches[edge];
            for (int i = 0; i < num_channels; ++i) { 
                add_val(vals[i], i, starting_pos, feature_caches);
            }           
        } else {
            std::vector<void*>& feature_caches = create_cache(edge);
            for (int i = 0; i < num_channels; ++i) { 
                add_val(vals[i], i, starting_pos, feature_caches);
            }        
        } 
    }

    void mv_features(RagEdge<Label>* edge2, RagEdge<Label>* edge1)
    {
        edge1->set_size(edge2->get_size());
        edge_caches[edge1] = edge_caches[edge2];
        edge_caches.erase(edge2);        
    } 

    void merge_features(RagNode<Label>* node1, RagNode<Label>* node2);
    void merge_features2(RagNode<Label>* node1, RagNode<Label>* node2, RagEdge<Label>* edge );
    void merge_features(RagEdge<Label>* edge1, RagEdge<Label>* edge2);


    void set_classifier(EdgeClassifier* pclfr){
    
     	eclfr = pclfr;
	std::vector<unsigned int> ignore_list;
	eclfr->get_ignore_featlist(ignore_list);
	for(size_t i=0; i<ignore_list.size(); i++)
	    ignore_set.insert(ignore_list[i]);
	
	printf("ignore features:");
	std::set<unsigned int>::iterator iit = ignore_set.begin();
	for(; iit != ignore_set.end(); iit++)
	    printf("%u ", *iit);
	printf("\n");
    }


    EdgeClassifier* get_classifier(){
	return eclfr;
    }  	

    //void set_tree_weights(std::vector<double>& pwts){
    //	eclfr->set_tree_weights(pwts);
    //}

    void add_median_feature();

#ifdef SETPYTHON
    void add_hist_feature(unsigned int num_bins, boost::python::list percentiles, bool use_diff);
#else
    void add_hist_feature(unsigned int num_bins, std::vector<double> percentiles, bool use_diff);
#endif
    void add_moment_feature(unsigned int num_moments, bool use_diff);
    
    void add_inclusiveness_feature(bool use_diff);

    double get_prob(RagEdge<Label>* edge);
    std::vector<double> get_prob_batch(std::vector< RagEdge<Label>*> edge);

    void get_responses(RagEdge<Label>* edge, vector<double>& responses);

    void compute_all_features(RagEdge<Label>* edge, std::vector<double>&);
    void compute_node_features(RagNode<Label>* edge, std::vector<double>&);
    
    void find_useless_features(std::vector< std::vector<double> >& all_features);

    ~FeatureMgr();


    void copy_channel_features(FeatureMgr *pfmgr);   	

    void copy_cache(std::vector<void*>& src_edge_cache, RagEdge<Label>* edge);	
    void copy_cache(std::vector<void*>& src_node_cache, RagNode<Label>* node);	

    void print_cache(RagEdge<Label>* edge);	
    void print_cache(RagNode<Label>* node);	

    std::vector<std::vector<FeatureCompute*> >& get_channel_features(){return channels_features;}; 	
    EdgeCaches& get_edge_cache(){return edge_caches;};  		
    NodeCaches& get_node_cache(){return node_caches;};  		
  private:
    void compute_diff_features(std::vector<void*>* caches1, std::vector<void*>* caches2, std::vector<double>& feature_results, RagEdge<Label>* edge);
    void compute_diff_features2(std::vector<void*>* caches1, std::vector<void*>* caches2, std::vector<double>& feature_results, RagEdge<Label>* edge);
     
    void compute_features(unsigned int prediction_type, std::vector<void*>* caches, std::vector<double>& feature_results, RagEdge<Label>* edge, unsigned int node_num);
    void compute_features2(unsigned int prediction_type, std::vector<void*>* caches, std::vector<double>& feature_results, RagEdge<Label>* edge, unsigned int node_num);

    void add_val(unsigned char val, unsigned int channel, unsigned int& starting_pos, std::vector<void *>& feature_caches)
    {
        std::vector<FeatureCompute*>& features = channels_features[channel];
        for (int i = 0; i < features.size(); ++i) {
            features[i]->add_point(val, feature_caches[starting_pos], channel > 0); 
            ++starting_pos;
        }
    }
   
    // !! assume all edge/node caches
    std::vector<void*>& create_cache(RagEdge<Label>* edge)
    {
        edge_caches[edge] = std::vector<void*>();
        std::vector<void*>& caches = edge_caches[edge];
        unsigned int pos = 0;
        for (unsigned int i = 0; i < num_channels; ++i) {
            std::vector<FeatureCompute*>& features = channels_features[i];
            for (int j = 0; j < features.size(); ++j) {
                caches.push_back(features[j]->create_cache());
                ++pos;
            } 
        }
        return caches;
    }
    std::vector<void*>& create_cache(RagNode<Label>* node)
    {
        node_caches[node] = std::vector<void*>();
        std::vector<void*>& caches = node_caches[node];
        unsigned int pos = 0;
        for (unsigned int i = 0; i < num_channels; ++i) {
            std::vector<FeatureCompute*>& features = channels_features[i];
            for (int j = 0; j < features.size(); ++j) {
                caches.push_back(features[j]->create_cache());
                ++pos;
            } 
        }
        return caches;
    }    

    void add_feature(unsigned int channel, FeatureCompute * feature, std::vector<bool>& feature_modes);


    EdgeCaches edge_caches;
    NodeCaches node_caches;
    unsigned int num_features;

    bool specified_features;
    unsigned int num_channels;
    bool simulate;

    std::vector<std::vector<FeatureCompute*> > channels_features;
    std::vector<std::vector<FeatureCompute*> > channels_features_equal;
    std::vector<std::vector<std::vector<bool> > > channels_features_modes;

#ifdef SETPYTHON
    boost::python::object pyfunc;
#endif
    bool has_pyfunc;


    EdgeClassifier* eclfr;	 
    std::set<unsigned int> ignore_set;
};




}

#endif


