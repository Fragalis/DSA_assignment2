#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <list>

using namespace std;

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
};

string printkDTreeNode(const kDTreeNode &node) {
    stringstream ss;
    ss << "(";
    if(node.data.size() > 0) {
        ss << node.data[0];
        for(int i = 1; i < node.data.size(); ++i) {
            ss << ",";
            ss << node.data[i];
        }
    }
    ss << ")";
    return ss.str();
}

void print_tree(kDTreeNode *node, string prefix) {
    if(!node) {
        cout << "EMPTY" << endl;
        return;
    }
    if(prefix == "") cout << "ROOT: ";
    else cout << prefix << " ";
    cout << printkDTreeNode(*node);
    cout << "\n";
    if(node->left) print_tree(node->left, prefix + "0");
    if(node->right) print_tree(node->right, prefix + "1");
}

class kDTree
{
private:
    int k;
    // Helpers:
    void clear() const;
    void clear_helper(kDTreeNode *node) const;
    void copy_helper(kDTreeNode *node);
    
    string inorder_helper(kDTreeNode *node) const;
    string preorder_helper(kDTreeNode *node) const;
    string postorder_helper(kDTreeNode *node) const;

    int height_record(kDTreeNode *node) const;
    int count_helper(kDTreeNode *node) const;
    int leaf_count_helper(kDTreeNode *node) const;

    kDTreeNode *insert_helper(kDTreeNode *node, const vector<int> &point, int idx);
    kDTreeNode *find_replacement_helper(kDTreeNode *node, int divisor, int idx, int &dim);
    kDTreeNode *find_parent_node_helper(kDTreeNode *node, const kDTreeNode &finder, int &dir);
    kDTreeNode *delete_helper(kDTreeNode *node, const vector<int> &point, int idx);
    bool search_helper(kDTreeNode *node, const vector<int> &point, int idx);

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

kDTree::kDTree(int k)
{
    this->k = k;
    this->root = nullptr;
};

void kDTree::clear() const {
    clear_helper(this->root);
};
void kDTree::clear_helper(kDTreeNode *node) const {
    if(!node) return;
    clear_helper(node->left);
    node->left = NULL;
    clear_helper(node->right);
    node->right = NULL;
    delete node;
};
kDTree::~kDTree()
{
    this->k = 0;
    this->clear();
    this->root = NULL;
};

const kDTree &kDTree::operator=(const kDTree &other)
{
    if(this != &other) {
        kDTree temp(other);
        swap(temp.k, this->k);
        swap(temp.root, this->root);
    }
    return *this;
};

void kDTree::copy_helper(kDTreeNode *node) {
    this->insert(node->data);
    if(node->left) copy_helper(node->left);
    if(node->right) copy_helper(node->right);
};
kDTree::kDTree(const kDTree &other)
{
    this->k = other.k;
    this->root = NULL;
    kDTreeNode *traverse = other.root;
    copy_helper(traverse);
};

string kDTree::inorder_helper(kDTreeNode *node) const {
    if(node == NULL) return "";
    string res = "";
    if(node->left) res += (inorder_helper(node->left) + " ");
    res += printkDTreeNode(*node);
    if(node->right) res += (" " + inorder_helper(node->right));
    return res;
}
void kDTree::inorderTraversal() const
{
    cout << inorder_helper(this->root);
};

string kDTree::preorder_helper(kDTreeNode *node) const {
    if(node == NULL) return "";
    string res = "";
    res += printkDTreeNode(*node);
    if(node->left) res += (" " + preorder_helper(node->left));
    if(node->right) res += (" " + preorder_helper(node->right));
    return res;
};
void kDTree::preorderTraversal() const
{
    cout << preorder_helper(this->root);
};

string kDTree::postorder_helper(kDTreeNode *node) const {
    if(node == NULL) return "";
    string res = "";
    if(node->left) res += (postorder_helper(node->left) + " ");
    if(node->right) res += (postorder_helper(node->right) + " ");
    res += printkDTreeNode(*node);
    return res;
};
void kDTree::postorderTraversal() const
{
    cout << postorder_helper(this->root);
};

int kDTree::height_record(kDTreeNode *node) const {
    if(node == NULL) return 0;
    return 1 + max(height_record(node->left), height_record(node->right));
};
int kDTree::height() const
{
    return height_record(this->root);
};

int kDTree::count_helper(kDTreeNode *node) const {
    if(node == NULL) return 0;
    return 1 + count_helper(node->left) + count_helper(node->right);
};
int kDTree::nodeCount() const
{   
    return count_helper(this->root);
};

int kDTree::leaf_count_helper(kDTreeNode *node) const {
    if(node == NULL) return 0;
    if(node->left == NULL && node->right == NULL) return 1;
    return leaf_count_helper(node->left) + leaf_count_helper(node->right);
};
int kDTree::leafCount() const
{
    return leaf_count_helper(this->root);
};


kDTreeNode *kDTree::insert_helper(kDTreeNode *node, const vector<int> &point, int idx) {
    if(node == NULL) {
        kDTreeNode *new_node = new kDTreeNode(point);
        return new_node;
    }
    // insert left
    if(point[idx] < node->data[idx]) node->left = insert_helper(node->left, point, (idx + 1) % (this->k));
    // insert right
    else node->right = insert_helper(node->right, point, (idx + 1) % (this->k));
    return node;
};
void kDTree::insert(const vector<int> &point)
{
    // if the point's length is not k -> ignore
    if(this->k != point.size()) return;
    this->root = insert_helper(this->root, point, 0);
};

kDTreeNode *kDTree::find_parent_node_helper(kDTreeNode *node, const kDTreeNode &finder, int &dir) {
    if(!node) return NULL;
    if(node->left && node->left == &finder) {
        dir = -1;
        return node;
    }
    if(node->right && node->right == &finder) {
        dir = 1;
        return node;
    }
    kDTreeNode *left = find_parent_node_helper(node->left, finder, dir);
    kDTreeNode *right = find_parent_node_helper(node->right, finder, dir);
    return (left)? left : right;
};
kDTreeNode *kDTree::find_replacement_helper(kDTreeNode *node, int divisor, int idx, int &dim) {
    if(idx == divisor) {
        if(node->left == NULL) return node;
        int next_dim = (dim + 1) % (this->k);
        return find_replacement_helper(node->left, divisor, (idx + 1) % (this->k), next_dim);
    }
    int r_dim = (dim + 1) % (this->k);
    int l_dim = (dim + 1) % (this->k);
    kDTreeNode *left_node = NULL;
    kDTreeNode *right_node = NULL;
    if(node->left) left_node = find_replacement_helper(node->left, divisor, (idx + 1) % (this->k), l_dim);
    if(node->right) right_node = find_replacement_helper(node->right, divisor, (idx + 1) % (this->k), r_dim);

    int dim_curr = node->data[divisor];
    int dim_left = (left_node)? left_node->data[divisor] : INT_MAX;
    int dim_right = (right_node)? right_node->data[divisor] : INT_MAX;

    if(dim_right < dim_left && dim_right < dim_curr) {
        dim = r_dim;
        return right_node;
    }
    if(dim_left <= dim_right && dim_left < dim_curr) {
        dim = l_dim;
        return left_node;
    }
    return node;
};
kDTreeNode *kDTree::delete_helper(kDTreeNode *node, const vector<int> &point, int idx) {
    if(node == NULL) return NULL;
    // if on the left
    if(point[idx] < node->data[idx]) {
        node->left = delete_helper(node->left, point, (idx + 1) % (this->k));
        return node;
    }
    if(point[idx] > node->data[idx]) {
        node->right = delete_helper(node->right, point, (idx + 1) % (this->k));
        return node;
    }
    // assume we found the node
    bool find_flag = true;
    // if not at left
    if(point[idx] == node->data[idx]) {
        for(unsigned int i = 0; i < this->k; ++i) {
            if(point[i] != node->data[i]) {
                // if it's not the node
                find_flag = false;
                break;
            }
        }
    }
    if(!find_flag) node->right = delete_helper(node->right, point, (idx + 1) % (this->k));
    else {
        // leaf case
        if(node->left == NULL && node->right == NULL) {
            delete node;
            return NULL;
        }

        // right child case
        if(node->right != NULL) {
            int divisor = idx;
            int dim = (idx + 1) % (this->k);
            kDTreeNode *replacement = find_replacement_helper(node->right, divisor, (idx + 1) % (this->k), dim);
            for(int i = 0; i < this->k; ++i) node->data[i] = replacement->data[i];
            int dir = 0;
            kDTreeNode *pre_replacement = find_parent_node_helper(node, *replacement, dir);
            if(dir == -1) {
                pre_replacement->left = delete_helper(replacement, node->data, dim);
            }
            if(dir == 1) {
                pre_replacement->right = delete_helper(replacement, node->data, dim);
            }
        }
        // only left child case
        else {
            int divisor = idx;
            int dim = (idx + 1) % (this->k);
            kDTreeNode *replacement = find_replacement_helper(node->left, divisor, (idx + 1) % (this->k), dim);
            for(int i = 0; i < this->k; ++i) node->data[i] = replacement->data[i];
            swap(node->right, node->left);
            int dir = 0;
            kDTreeNode *pre_replacement = find_parent_node_helper(node, *replacement, dir);
            if(dir == -1) {
                pre_replacement->left = delete_helper(replacement, node->data, dim);
            }
            if(dir == 1) {
                pre_replacement->right = delete_helper(replacement, node->data, dim);
            }
        }
    }
    return node;
};
void kDTree::remove(const vector<int> &point)
{
    if(this->k != point.size()) return;
    this->root = delete_helper(this->root, point, 0);
};

bool kDTree::search_helper(kDTreeNode *node, const vector<int> &point, int idx) {
    if(node == NULL) return false;
    // if on the left
    if(point[idx] < node->data[idx]) return search_helper(node->left, point, (idx + 1) % (this->k));
    // assume we found the node
    bool find_flag = true;
    // if not at left
    if(point[idx] == node->data[idx]) {
        for(unsigned int i = 0; i < this->k; ++i) {
            if(point[i] != node->data[i]) {
                // if it's not the node
                find_flag = false;
                break;
            }
        }
    }
    // it can be this node or the node the node on the right
    return find_flag || search_helper(node->right, point, (idx + 1) % (this->k));
};
bool kDTree::search(const vector<int> &point)
{
    if(this->k != point.size()) return false;
    return search_helper(this->root, point, 0);
};

void kDTree::buildTree(const vector<vector<int>> &pointList)
{
    return;
};
void kDTree::neighbour_finder(kDTreeNode *node, const vector<int> &target, kDTreeNode *&best, int idx, long &best_distance) {
    if(!node) return;
    // cout << "CURR : " << printkDTreeNode(*node) << " BEST : " << printkDTreeNode(*best) << endl;
    int dir = 0;
    if(target[idx] < node->data[idx]) {
        if(!node->left) {
            for(int i = 0; i < this->k; ++i)
                best_distance += (long)abs(pow(target[i] - node->data[i], 2));
            best = node;
            return;
        }
        neighbour_finder(node->left, target, best, (idx + 1) % (this->k), best_distance);
        dir = 1;
    }
    else {
        if(!node->right) {
            for(int i = 0; i < this->k; ++i)
                best_distance += (long)abs(pow(target[i] - node->data[i], 2));
            best = node;
            return;
        }
        neighbour_finder(node->right, target, best, (idx + 1) % (this->k), best_distance);
        dir = -1;
    }

    long curr_distance = 0;
    for(int i = 0; i < this->k; ++i)
        curr_distance += (long)abs(pow(target[i] - node->data[i], 2));
    if(curr_distance < best_distance) {
        best_distance = curr_distance;
        best = node;
    }
    long plane_distance = abs(pow(target[idx] - node->data[idx], 2));
    if(best_distance >= plane_distance) {
        best_distance = plane_distance;
        if(dir == -1) neighbour_finder(node->left, target, best, (idx + 1) % (this->k), best_distance);
        if(dir == 1) neighbour_finder(node->right, target, best, (idx + 1) % (this->k), best_distance);
    }
};
void kDTree::nearestNeighbour(const vector<int> &target, kDTreeNode *&best)
{
    long best_distance = 0;
    best = this->root;
    neighbour_finder(this->root, target, best, 0, best_distance);
    // cout << printkDTreeNode(*best) << endl;
};

void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList)
{
    kDTree temp(*this);
    int number_of_neighbour = 0;
    kDTreeNode *neighbour = temp.root;
    while(number_of_neighbour < k) {
        temp.nearestNeighbour(target, neighbour);
        bestList.push_back(neighbour);
        temp.remove(neighbour->data);
        neighbour = temp.root;
        ++number_of_neighbour;
    }
};

void test_print_node() {
    vector<int> test = {};
    kDTreeNode t(test);
    cout << printkDTreeNode(t);
}

void test_solo_tree() {
    vector< vector<int> > store = 
    {
        {5, 6}, {2, 2}, {7, 3}, {2, 8}, {8, 1}, {9, 2}, {3, 5}
    };
    int k = store[0].size();
    kDTree *tree = new kDTree(k);
    for(auto &v : store) {
        tree->insert(v);
    }

    print_tree(tree->root, "");
    // tree->preorderTraversal();
    cout << endl;
    // tree->inorderTraversal();
    cout << endl;
    // tree->postorderTraversal();
    cout << endl;
    for(auto &v : store) {
        tree->remove(v);
        print_tree(tree->root, "");
    }
}

void test_copy_tree() {
    vector< vector<int> > store = 
    {
        {5, 6}, {2, 2}, {7, 3}, {2, 8}, {8, 1}, {9, 2}, {3, 5}
    };
    int k = store[0].size();
    kDTree *tree = new kDTree(k);
    for(auto &v : store) {
        tree->insert(v);
    }
    cout << "Original: ";
    tree->inorderTraversal();
    cout << endl;

    // copy constructor
    cout << "Copy Constructor: ";
    kDTree copy_tree(*tree);
    copy_tree.inorderTraversal();
    cout << endl;

    // assignment operator
    cout << "Assignment Operator: ";
    kDTree ass_tree;
    ass_tree = *tree;
    ass_tree.inorderTraversal();
    cout << endl;

    // destructor
    cout << "Destructor: ";
    tree->~kDTree();
    tree->inorderTraversal();
    cout << endl;

    cout << "Copy Constructor (expect an original tree here): ";
    copy_tree.inorderTraversal();
    cout << endl;

    cout << "Assignment Operator (expect an original tree here): ";
    ass_tree.inorderTraversal();
    cout << endl;
}

void test_methods() {
    vector< vector<int> > store = 
    {
        {5, 6}, {2, 2}, {7, 3}, {2, 8}, {8, 1}, {9, 2}, {3, 5}
    };
    int k = store[0].size();
    kDTree *tree = new kDTree(k);
    for(auto &v : store) {
        tree->insert(v);
    }
    cout << "Original: ";
    tree->inorderTraversal();
    cout << endl;

    // 4
    cout << tree->height() << endl;
    // 7
    cout << tree->nodeCount() << endl;
    // 2
    cout << tree->leafCount() << endl;
}

void test_build_tree() {
    return;
}

void test_neighbours() {
    vector< vector<int> > store = 
    {
        {5, 6}, {2, 2}, {7, 3}, {2, 8}, {8, 2}, {8, 7}, {9, 2}, {3, 5}, {9, 5}
    };
    int k = store[0].size();
    kDTree *tree = new kDTree(k);
    for(auto &v : store) {
        tree->insert(v);
    }
    print_tree(tree->root, "");
    vector<int> target = {9,3};
    kDTreeNode *neighbour;
    tree->nearestNeighbour(target, neighbour);
    cout << "nearestNeighbour TEST: " << printkDTreeNode(*neighbour) << endl;

    vector<kDTreeNode*> neighbour_list;
    tree->kNearestNeighbour(target, 3, neighbour_list);
    cout << "kNearestNeighbour TEST: " << endl;
    for(kDTreeNode *node : neighbour_list) {
        cout << printkDTreeNode(*node);
        cout << endl;
    }
}

int main() {
    // test_print_node();
    // test_solo_tree();
    // test_copy_tree();
    // test_methods();
    // test_build_tree();
    test_neighbours();
}