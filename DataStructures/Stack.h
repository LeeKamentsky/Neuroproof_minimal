#include <vector>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include "Rag.h"
#include "AffinityPair.h"
#include "../FeatureManager/FeatureManager.h"
#include "Glb.h"
// #include <boost/python.hpp>
// #include <boost/tuple/tuple.hpp>
#include "../Algorithms/MergePriorityFunction.h"
#include "../Algorithms/MergePriorityQueue.h"
#include <time.h>

//#include "../Algorithms/BatchMergeMRFh.h"
#include "../Utilities/unique_row_matrix.h"

#include "../Watershed/vigra_watershed.h"

#include "../Algorithms/RagAlgs.h"


#ifndef STACK_H
#define STACK_H

namespace NeuroProof {
class LabelCount;
class Stack {
  public:
    Stack(Label* watershed_, int depth_, int height_, int width_, int padding_=1) : watershed(watershed_), depth(depth_), height(height_), width(width_), padding(padding_), feature_mgr(0), median_mode(false), prediction_array(0), gtruth(NULL), merge_mito(true)
    {
        rag = new Rag<Label>();
    }

    void build_rag();

    void merge_rags (Rag<Label>* rag1, Rag<Label>* rag2, FeatureMgr* fm1, FeatureMgr* fm2);

    void build_rag_recurse (Rag<Label>* rag1, FeatureMgr* fm1,  
                                                    int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, bool use_mito_prob);

    Label * get_label_volume();

    void determine_edge_locations();

//     bool add_edge_constraint(boost::python::tuple loc1, boost::python::tuple loc2);

    Label get_body_id(unsigned int x, unsigned int y, unsigned int z);


    //boost::python::tuple get_edge_loc(RagEdge<Label>* edge);
//     void get_edge_loc(RagEdge<Label>* edge, Label& x, Label& y, Label& z);



    void add_prediction_channel(unsigned char* prediction_array_)
    {
        prediction_array.push_back(prediction_array_);
    
        if (feature_mgr) {
            feature_mgr->add_channel();
        }
    }

    bool is_orphan(RagNode<Label>* node)
    {
        return !(node->is_border());
    }
        
    Rag<Label>* get_rag()
    {
        return rag;
    }


    void build_rag_loop(Rag<Label>* rag, FeatureMgr* feature_man, 
                                            int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, bool use_mito_prob);
    void build_rag_loop_tfk(Rag<Label>* rag, FeatureMgr* feature_man, 
					    int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, bool use_mito_prob);
    void build_rag_loop_leek(Rag<Label>* rag, FeatureMgr* feature_man, 
					    int x_start, int x_end, int y_start, int y_end, int z_start, int z_end, bool use_mito_prob);
    double get_edge_weight(RagEdge<Label>* edge)
    {
        return feature_mgr->get_prob(edge);
    }

    void set_basic_features();
	

    FeatureMgr * get_feature_mgr()
    {
        return feature_mgr;
    }
    void set_feature_mgr(FeatureMgr * feature_mgr_)
    {
        feature_mgr = feature_mgr_;
    }


    int get_num_bodies()
    {
        return rag->get_num_regions();        
    }

    int get_width() const
    {
        return width;
    }

    int get_height() const
    {
        return height;
    }

    int get_depth() const
    {
        return depth;
    }

    int remove_inclusions();

    ~Stack()
    {
        delete rag;
        if (!prediction_array.empty()) {
            for (unsigned int i = 0; i < prediction_array.size(); ++i) {
                delete prediction_array[i];
            }
        }

        if (feature_mgr) {
            delete feature_mgr;
        }
    }
    struct DFSStack {
        Label previous;
        RagNode<Label>* rag_node;  
        int count;
        int start_pos;
    };

    void set_groundtruth(Label* pgt) {gtruth = pgt; }
    void compute_groundtruth_assignment();     			
    void compute_contingency_table();
    void compute_vi();
    void modify_assignment_after_merge(Label node_keep, Label node_remove);
//     void write_graph(string);
    int decide_edge_label(RagNode<Label>* node1, RagNode<Label>* node2);
    void set_merge_mito(bool flag, double pthd=0.35)
    {
	merge_mito=flag;
	mito_thd = pthd;
    };

  protected:
    void biconnected_dfs(std::vector<DFSStack>& dfs_stack);
    
    Rag<Label> * rag;
    Label* watershed;
    std::vector<unsigned char*> prediction_array;
    std::tr1::unordered_map<Label, Label> watershed_to_body; 
    std::tr1::unordered_map<Label, std::vector<Label> > merge_history; 
    //typedef std::tuple<unsigned int, unsigned int, unsigned int> Location;
//     typedef boost::tuple<unsigned int, unsigned int, unsigned int> Location;
    typedef std::tr1::unordered_map<RagEdge<Label>*, unsigned long long> EdgeCount; 
//     typedef std::tr1::unordered_map<RagEdge<Label>*, Location> EdgeLoc; 
    
    int depth, height, width;
    int padding;

//     boost::shared_ptr<PropertyList<Label> > node_properties_holder;
    std::tr1::unordered_set<Label> visited;
    std::tr1::unordered_map<Label, int> node_depth;
    std::tr1::unordered_map<Label, int> low_count;
    std::tr1::unordered_map<Label, Label> prev_id;
    std::vector<std::vector<OrderedPair<Label> > > biconnected_components; 
    
    std::vector<OrderedPair<Label> > stack;
   
    EdgeCount best_edge_z;
//     EdgeLoc best_edge_loc;

    FeatureMgr * feature_mgr;
    bool median_mode;

    bool merge_mito;
    double mito_thd;

    Label* gtruth; 	// both derived classes may use groundtruth for either learning or validation
    std::multimap<Label, Label>	assignment;
    std::multimap<Label, std::vector<LabelCount> > contingency;	



};



class StackLearn: public Stack{

public:

    StackLearn(Label* watershed_, int depth_, int height_, int width_, int padding_=1): Stack(watershed_, depth_, height_, width_, padding_){}

    void learn_edge_classifier_lash(double, UniqueRowFeature_Label& , std::vector<int>&, bool prune_feature, const char* clfr_filename = NULL);
    void learn_edge_classifier_queue(double, UniqueRowFeature_Label& , std::vector<int>&, bool prune_feature=false, const char* clfr_filename = NULL);
    void learn_edge_classifier_flat(double, UniqueRowFeature_Label& , std::vector<int>&, bool prune_feature=false, const char* clfr_filename = NULL);
    void learn_edge_classifier_flat_subset(double, UniqueRowFeature_Label& , std::vector<int>&, const char* clfr_filename = NULL);

};

class StackPredict: public Stack{

public:

    StackPredict(Label* watershed_, int depth_, int height_, int width_, int padding_=1): Stack(watershed_, depth_, height_, width_, padding_){}
    
    void agglomerate_rag(double threshold, bool use_edge_weight=false, string output_path="", string classifier_path="");
    void agglomerate_rag_queue(double threshold, bool use_edge_weight=false, string output_path="", string classifier_path="");     			
    void agglomerate_rag_flat(double threshold, bool use_edge_weight=false, string output_path="", string classifier_path="");
    void agglomerate_rag_mrf(double threshold, bool read_off, string output_path, string classifier_path);
    void agglomerate_rag_size(double threshold);

    void merge_mitochondria_a();
    void absorb_small_regions(float* prediction_vol, Label* label_vol);
    void absorb_small_regions2(float* prediction_vol, Label* label_vol, size_t);

};

class LabelCount{
public:
      Label lbl;
      size_t count;
      LabelCount(): lbl(0), count(0) {};	 	 		
      LabelCount(Label plbl, size_t pcount): lbl(plbl), count(pcount) {};	 	 		
};


}

#endif
