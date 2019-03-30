/**
 * HEADER todo
 */

#include <iostream>
#include <vector>

using namespace std;

// selectionSort
//
// Effectue le tri par sélection des éléments entre begin
// et end (non inclus). Doit appeler display() après chaque
// échange.
//
// A COMPLETER
template < typename RandomAccessIterator >
void selectionSort(RandomAccessIterator begin, RandomAccessIterator end) {

    for (RandomAccessIterator i = begin; i < end - 1; ++i) {
        RandomAccessIterator imin = i;

        for (RandomAccessIterator j = i + 1; j < end; ++j) {

            if (*j < *imin) {
                imin = j;
            }
        }

        swap(*i, *imin);
    }
}


// selectPivot(begin,end)
//
// choisit un pivot pour le tri rapide parmi le tableau
// entre begin et end (non inclus). Calcule la médiane
// entre le premier, le dernier et l'élément central.
// retourne un iterateur du même type que begin et end
// pointant sur la valeur choisie.
//
// NE RIEN MODIFIER DANS CETTE FONCTION

template < typename RandomAccessIterator >
RandomAccessIterator selectPivot(const RandomAccessIterator begin, const RandomAccessIterator end) {
    RandomAccessIterator p1 = begin;
    RandomAccessIterator p2 = begin + ( end-begin ) / 2;
    RandomAccessIterator p3 = end-1;

    if(*p1 < *p2) {
        if( *p2 < *p3 ) return p2;
        else if(*p1 < *p3) return p3;
        return p1;
    } else {
        if( *p1 < *p3 ) return p1;
        else if(*p2 < *p3) return p3;
        return p2;
    }
}

/**
 * @brief We create a partition, all elements lower than the pivot will be place on the left of it and all
 *        elements greather than the pivot will be place on the right of it.
 *
 * @tparam RandomAccessIterator
 * @param begin The begin iterator
 * @param end The end iterator that point after the last element
 * @return An iterator on the pivot's position in the partition. This element is the only one for which we are sure
 *         that it have the right position.
 */
template < typename RandomAccessIterator >
RandomAccessIterator partition(RandomAccessIterator begin, RandomAccessIterator end) {
    RandomAccessIterator pivot = end - 1; // -1 because the end iterator point after the last element
    RandomAccessIterator i = begin - 1; // -1 because we use do while, so he will be directly increments before the loop condition
    RandomAccessIterator j = pivot;

    while(true) {
        // Increment i while we can
        do {
            ++i;
        } while (*i < *pivot);

        // Increment j while we can
        do {
            --j;
        } while (j > begin and *pivot < *j);

        // stop the loop when iterators meet each other
        if (i >= j) {
            break;
        }

        swap(*i, *j);
    }

    swap(*i, *pivot);

    return i;
}

// quickSort
//
// Effectue le tri rapide des éléments entre begin
// et end (non inclus). Doit appeler selectPivot(...)
// pour le choix du pivot, et display() après chaque
// partition
//
// A COMPLETER

template < typename RandomAccessIterator >
void quickSort(RandomAccessIterator begin, RandomAccessIterator end) {
    // trivial case
    if(end - 1 <= begin) {
        return;
    }

    // Pivot selection
    RandomAccessIterator p = selectPivot(begin, end);
    // Swap the last element with the pivot before the partition call
    swap(*(end-1), *p);
    RandomAccessIterator i = partition(begin, end);

    // Our two recursive calls, one for the left part, one for the right part. The element in the middle
    // (return by the partition) is correctly placed
    quickSort(begin, i);
    quickSort(i+1, end);
}

template < typename RandomAccessIterator >
void display(RandomAccessIterator begin, RandomAccessIterator end) {
    cout << "[";
    for(auto i = begin; i != end; ++i) {
        if (i != begin){
            cout << ", ";
        }
        cout << *i;
    }
    cout << "]" << endl;
}



int main(int argc, const char * argv[]) {

    vector<int> vect {7,4,5,2,9,8};

    display(vect.begin(), vect.end());
    selectionSort(vect.begin(), vect.end());
    display(vect.begin(), vect.end());

    return 0;
}