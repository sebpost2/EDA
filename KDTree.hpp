// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <size_t N, typename ElemType>
class KDTree {
 public:

     struct node {
         node* nodes[2];
         int height;
         ElemType value;
         Point<N> point;
         node() {
             nodes[0] = nullptr;
             nodes[1] = nullptr;
             point(0, 0);
             value = 0;
             height = 0;
         }
         node(const Point<N> val) {
             nodes[0] = nullptr;
             nodes[1] = nullptr;
             point=val;
             height = 0;
             value = 0;
         }

         node(const Point<N> val, const int height) {
             nodes[0] = nullptr;
             nodes[1] = nullptr;
             point = val;
             this->height = height;
             value = 0;
         }

         node(const Point<N> val, const int height, const ElemType value) {
             nodes[0] = nullptr;
             nodes[1] = nullptr;
             point = val;
             this->height = height;
             this->value = value;
         }
     };

  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);

  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  node* find(node* it, const Point<N>& val) const;

  node* dfsImitate(node* it);

  void destroy(node* it);

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

 private:
     node* root;
     size_t dimension_;
     size_t size_;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    root = NULL;
    size_=0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    destroy(root);
}

template <std::size_t N, typename ElemType>
void KDTree<N, ElemType>::destroy( node* it)  {
    if (it == NULL) {
        return;
    }
    if (it->nodes[0] != NULL) {
        destroy(it->nodes[0]);
    }
    if (it->nodes[1] != NULL) {
        destroy(it->nodes[1]);
    }
    delete it;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
    root = rhs.root;
    size_ = rhs.size_;
}

template <std::size_t N, typename ElemType>
typename KDTree<N, ElemType>::node* KDTree<N, ElemType>::find(node* it, const Point<N>& pt) const  {
    if (it == NULL || it->point == pt) 
        return it;

    if (pt[it->height % N] < it->point[it->height % N]) {
        if (it->nodes[0] == NULL)
            return it; 
        return find(it->nodes[0], pt);
    }
    else {
        if (it->nodes[1] == NULL)
            return it;
        return find(it->nodes[1], pt);
    }
}

template <std::size_t N, typename ElemType>
typename KDTree<N, ElemType>::node* KDTree<N, ElemType>::dfsImitate(typename KDTree<N, ElemType>::node* it) {
    node* temp = new node;
    temp = it;
    temp->nodes[0] = dfsImitate(it->nodes[0]);
    temp->nodes[1] = dfsImitate(it->nodes[1]);
    return temp;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
    if (root != rhs.root) {
        root = rhs.root;
        size_ = rhs.size_;
  }
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    if (root == NULL) {
        return true;
  }
  return false;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
    node* Temp = find(root,pt);
    if (Temp->point == pt) {
        return true;
    }
    return false;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    node* Temp = find(root,pt);

    if (Temp == NULL) {
        root = new node(pt, 0, value);
        size_++;
        return;
    }

    if (Temp->point == pt) {
        return;
    }

    int height = Temp->height;
    size_++;
    node* son = new node(pt, height + 1, value);
    if (pt[height % N] < Temp->point[height % N]) {
        Temp->nodes[0] = son;
    }
    else {
        Temp->nodes[1] = son;
    }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    /*node* Temp = find(pt);
    if (Temp->point != pt || Temp == NULL) {
        throw std::out_of_range("Point not found in the KD-Tree");
    }
    else {
        return Temp->value;
    }*/
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
    const KDTree<N, ElemType>& CONST = *this;
    return const_cast<ElemType&>(CONST.at(pt));
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
    node* Temp = find(root,pt);
    if (Temp->point != pt || Temp == NULL) {
        throw std::out_of_range("Point not found in the KD-Tree");
    }
    else {
        return Temp->value;
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // TODO(me): Fill this in.
  ElemType new_element;
  return new_element;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  return values;
}

// TODO(me): finish the implementation of the rest of the KDTree class

#endif  // SRC_KDTREE_HPP_
