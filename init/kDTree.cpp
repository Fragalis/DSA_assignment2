#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */
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
    clear_helper(node->right);
    delete node;
    return;
};
kDTree::~kDTree()
{
    this->k = 0;
    clear();
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

void kDTree::inorder_helper(kDTreeNode *node) const {
    if(node == NULL) return;
    // string res = "";
    if(node->left) {
        inorder_helper(node->left);
        cout << " ";
    }
    cout << *node;
    if(node->right) {
        cout << " ";
        inorder_helper(node->right);
    }
}
void kDTree::inorderTraversal() const
{
    inorder_helper(this->root);
};

void kDTree::preorder_helper(kDTreeNode *node) const {
    if(node == NULL) return;
    // string res = "";
    cout << *node;
    if(node->left) {
        cout << " ";
        preorder_helper(node->left);
    }
    if(node->right) {
        cout << " ";
        preorder_helper(node->right);
    }
}
void kDTree::preorderTraversal() const
{
    preorder_helper(this->root);
};

void kDTree::postorder_helper(kDTreeNode *node) const {
    if(node == NULL) return;
    // string res = "";
    if(node->left) {
        postorder_helper(node->left);
        cout << " ";
    }
    if(node->right) {
        postorder_helper(node->right);
        cout << " ";
    }
    cout << *node;
}
void kDTree::postorderTraversal() const
{
    postorder_helper(this->root);
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
    if(point[idx] > node->data[idx]) return search_helper(node->right, point, (idx + 1) % (this->k));
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

void kDTree::merge(vector<vector<int>> &pointList, int left, int mid, int right, int idx) const {
    int const left_size = mid - left + 1;
    int const right_size = right - mid;
    vector<vector<int>> left_sublist;
    vector<vector<int>> right_sublist;

    for(int i = 0; i < left_size; ++i) left_sublist.push_back(pointList[left + i]);
    for(int i = 0; i < right_size; ++i) right_sublist.push_back(pointList[mid + 1 + i]);
    int left_idx = 0, right_idx = 0, merged_idx = left;
    while(left_idx < left_size && right_idx < right_size) {
        if(left_sublist[left_idx][idx] <= right_sublist[right_idx][idx]) {
            pointList[merged_idx] = left_sublist[left_idx];
            ++left_idx;
        }
        else {
            pointList[merged_idx] = right_sublist[right_idx];
            ++right_idx;
        }
        ++merged_idx;
    }
    while(left_idx < left_size) {
        pointList[merged_idx] = left_sublist[left_idx];
        ++left_idx;
        ++merged_idx;
    }
    while(right_idx < right_size) {
        pointList[merged_idx] = right_sublist[right_idx];
        ++right_idx;
        ++merged_idx;
    }
};
void kDTree::merge_sort_helper(vector<vector<int>> &pointList, int start, int end, int idx) const {
    if(start >= end) return;
    int mid = start + (end - start)/2;
    merge_sort_helper(pointList, start, mid, idx);
    merge_sort_helper(pointList, mid + 1, end, idx);
    // cout << "merge " << start << " to " << end << " with mid = " << mid << endl;
    merge(pointList, start, mid, end, idx);
};
kDTreeNode *kDTree::build_helper(vector<vector<int>> &pointList, int start, int end, int idx) {
    if(start > end) return NULL;
    // cout << "start = " << start << " end = " << end << endl;
    merge_sort_helper(pointList, start, end, idx);
    int mid = start + (end - start)/2;
    kDTreeNode *node = new kDTreeNode(pointList[mid]);
    node->left = build_helper(pointList, start, mid - 1, (idx + 1) % (this->k));
    node->right = build_helper(pointList, mid + 1, end, (idx + 1) % (this->k));
    return node;
};
void kDTree::buildTree(const vector<vector<int>> &pointList)
{
    vector<vector<int>> list = pointList;
    this->root = build_helper(list, 0, list.size() - 1, 0);
};

long kDTree::distance(const vector<int> &v1, const vector<int> &v2) {
    long best_distance = 0;
    for(int i = 0; i < this->k; ++i) best_distance += (long)abs(pow(v1[i] - v2[i], 2));
    return best_distance;
}

void kDTree::neighbour_finder(kDTreeNode *node, const vector<int> &target, kDTreeNode *&best, int idx) {
    if(!node) return;
    // cout << "CURR : " << printkDTreeNode(*node) << " BEST : " << printkDTreeNode(*best) << endl;
    int dir = 0;
    if(target[idx] < node->data[idx]) {
        if(!node->left) {
            if(!best || distance(target, node->data) < distance(target, best->data)) {
                best = node;
            }
            return;
        }
        neighbour_finder(node->left, target, best, (idx + 1) % (this->k));
        dir = 1;
    }
    else {
        if(!node->right) {
            if(!best || distance(target, node->data) < distance(target, best->data)) {
                best = node;
            }
            return;
        }
        neighbour_finder(node->right, target, best, (idx + 1) % (this->k));
        dir = -1;
    }

    long curr_distance = distance(target, node->data);
    long best_distance = distance(target, best->data);
    if(curr_distance < best_distance) {
        best = node;
    }

    long plane_distance = abs(pow(target[idx] - node->data[idx], 2));
    if(best_distance >= plane_distance) {
        if(dir != 1) neighbour_finder(node->left, target, best, (idx + 1) % (this->k));
        if(dir != -1) neighbour_finder(node->right, target, best, (idx + 1) % (this->k));
    }
};
void kDTree::nearestNeighbour(const vector<int> &target, kDTreeNode *&best)
{
    best = this->root;
    neighbour_finder(this->root, target, best, 0);
    // cout << printkDTreeNode(*best) << endl;
};

void kDTree::kNearestNeighbour(const vector<int> &target, int k, vector<kDTreeNode *> &bestList)
{
    kDTree temp(*this);
    int number_of_neighbour = 0;
    kDTreeNode *neighbour = temp.root;
    while(number_of_neighbour < k) {
        temp.nearestNeighbour(target, neighbour);
        kDTreeNode *return_val = new kDTreeNode(neighbour->data);
        bestList.push_back(return_val);
        temp.remove(neighbour->data);
        neighbour = temp.root;
        ++number_of_neighbour;
    }
};

kNN::kNN(int k) {
    this->k = k;
};
void kNN::fit(Dataset &X_train, Dataset &y_train) {
    this->X_train = &X_train;
    this->y_train = &y_train;
    vector<vector<int>> tree_data;
    for(auto lst = X_train.data.begin(); lst != X_train.data.end(); ++lst) {
        vector<int> temp;
        for(auto it = (*lst).begin(); it != (*lst).end(); ++it) temp.push_back(*it);
        tree_data.push_back(temp);
    }
    this->tree = new kDTree(tree_data[0].size());
    tree->buildTree(tree_data);
};
Dataset kNN::predict(Dataset &X_test) {
    Dataset y_pred;
    y_pred.columnName.push_back("label");

    for(auto lst = X_test.data.begin(); lst != X_test.data.end(); ++lst) {
        vector<int> target;
        for(auto it = (*lst).begin(); it != (*lst).end(); ++it) target.push_back(*it);
        vector<kDTreeNode*> best_neighbour;
        this->tree->kNearestNeighbour(target, this->k, best_neighbour);
        
        vector<int> count(10, 0);
        for(auto const &neighbour : best_neighbour) {
            int idx = 0;
            for(auto X_itr = this->X_train->data.begin(); X_itr != this->X_train->data.end(); ++X_itr) {
                vector<int> X_row;
                for(auto ele_itr = (*X_itr).begin(); ele_itr != (*X_itr).end(); ++ele_itr) {
                    X_row.push_back(*ele_itr);
                }
                if(X_row == neighbour->data) {
                    auto y_itr = this->y_train->data.begin();
                    advance(y_itr, idx);
                    count[*(*y_itr).begin()]++;
                }
                ++idx;
            }
        }

        int max_count = 0, num = 0;
        for(int pos = 0; pos < 10; ++pos) {
            // cout << count[pos] << " ";
            if(max_count < count[pos]) {
                num = pos;
                max_count = count[pos];
            }
        }
        // cout << endl;
        list<int> pred;
        pred.push_back(num);
        y_pred.data.push_back(pred);
    }
    return y_pred;
};
double kNN::score(const Dataset &y_test, const Dataset &y_pred) {
    int correct_count = 0;
    int test_count = 0;
    auto test_itr = y_test.data.begin();
    auto pred_itr = y_pred.data.begin();
    while(test_itr != y_test.data.end() && pred_itr != y_pred.data.end()) {
        ++test_count;
        if(*((*test_itr).begin()) == *((*pred_itr).begin())) ++correct_count;
        ++test_itr;
        ++pred_itr;
    }
    return (double)correct_count / (double)test_count;
};