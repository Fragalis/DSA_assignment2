#include "main.hpp"
#include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct kDTreeNode
{
    vector<int> data;
    kDTreeNode *left;
    kDTreeNode *right;
    kDTreeNode(vector<int> data, kDTreeNode *left = nullptr, kDTreeNode *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }
    
    friend ostream &operator<<(ostream &os, const kDTreeNode &node)
    {
        os << "(";
        for (int i = 0; i < node.data.size(); i++)
        {
            os << node.data[i];
            if (i != node.data.size() - 1)
            {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }
};

class kDTree
{
private:
    int k;

    // Constructor - Destructor Helpers
    void clear() const;
    void clear_helper(kDTreeNode *node) const;
    void copy_helper(kDTreeNode *node);
    
    // Printer Helpers
    void inorder_helper(kDTreeNode *node) const;
    void preorder_helper(kDTreeNode *node) const;
    void postorder_helper(kDTreeNode *node) const;

    // Misc Helpers
    int height_record(kDTreeNode *node) const;
    int count_helper(kDTreeNode *node) const;
    int leaf_count_helper(kDTreeNode *node) const;

    // Modify Helpers
    kDTreeNode *insert_helper(kDTreeNode *node, const vector<int> &point, int idx);
    kDTreeNode *find_parent_node_helper(kDTreeNode *node, const kDTreeNode &finder, int &dir);
    kDTreeNode *find_replacement_helper(kDTreeNode *node, int divisor, int idx, int &dim);
    kDTreeNode *delete_helper(kDTreeNode *node, const vector<int> &point, int idx);
    bool search_helper(kDTreeNode *node, const vector<int> &point, int idx);

    // Build Helpers
    void merge(vector<vector<int>> &pointList, int left, int mid, int right, int idx) const;
    void merge_sort_helper(vector<vector<int>> &pointList, int start, int end, int idx) const;
    kDTreeNode *build_helper(vector<vector<int>> &pointList, int start, int end, int idx);

    void neighbour_finder(kDTreeNode *node, const vector<int> &target, kDTreeNode *&best, int idx, long &best_distance);
public:
    kDTreeNode *root;
    kDTree(int k = 2);
    ~kDTree();

    const kDTree &operator=(const kDTree &other);
    kDTree(const kDTree &other);

    void inorderTraversal() const;
    void preorderTraversal() const;
    void postorderTraversal() const;
    int height() const;
    int nodeCount() const;
    int leafCount() const;

    void insert(const vector<int> &point);
    void remove(const vector<int> &point);
    bool search(const vector<int> &point);
    void buildTree(const vector<vector<int>> &pointList);
    void nearestNeighbour(const vector<int> &target, kDTreeNode *&best);
    void kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList);
};

class kNN
{
private:
    int k;
    Dataset *X_train;
    Dataset *y_train;
    int numClasses;
    kDTree *tree;

public:
    kNN(int k = 5);
    void fit(Dataset &X_train, Dataset &y_train);
    Dataset predict(Dataset &X_test);
    double score(const Dataset &y_test, const Dataset &y_pred);
};

// Please add more or modify as needed
