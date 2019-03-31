/**
 * HEADER todo
 */

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

enum class Sorts {
    SELECTION_SORT, QUICK_SORT, COUNTING_SORT
};

// functor for radix sort
template<typename T>
class Binary {
public:
    Binary(size_t shift) : shift(shift) {}


    /**
     * @TODO
     * @param value
     * @return
     */
    T operator()(T value) {
        unsigned mask = 0xff;

        return ((value >> 8 * shift) & mask);
    }

private:
    size_t shift;
};

/**
 * @brief Generic counting sort
 *
 * @note https://en.wikipedia.org/wiki/Counting_sort
 *
 * @param first    is the first element of the container to sort
 * @param last     is the last element of the container to sort
 * @param output   begining of the container where the result will be
 *                 written. Must be different from the original container
 * @param key  function that takes one element and return it's position within
 *             within the counting vector
 * @param max_key  max value that can be returned by key(). If it's value is -1
 *                 then the function needs to calculate it's value by looping through the
 *                 container to sort
 */
template<typename RandomAccessIterator, typename Fn>
void CountingSort(RandomAccessIterator first,
                  RandomAccessIterator last,
                  RandomAccessIterator output,
                  Fn key,
                  size_t max_key = -1) {


    // check if a max_key was given
    if (max_key == -1) {
        // if not, loop through the array to find if
        size_t max = 0;

        for (auto i = first; i != last; ++i) {

            if (key(*i) > max)
                max = key(*i);
        }

        max_key = max;
    }

    std::vector<size_t> count(max_key + 1, 0);

    // count the number of times each elements are found within the og container
    for (auto i = first; i != last; ++i) {
        auto val = key(*i);
        count[val] = count[val] + 1;
    }

    // loop through the cout vector to set the range in which elements will be placed̉
    for (size_t i = 1; i < count.size(); ++i) {
        count[i] = count[i] + count[i - 1];
    }

    for (auto i = last; i != first;) {
        // since we're looping from back to front, we need to
        //  decrement at the start of the loop and not at the end.
        // if we decrement at the end, there is a chance that random
        //  values get inserted into the ordered container
        --i;

        // get the real value of an element
        auto val = key(*i);
        // place it in the output in the correct area
        *(output + count[val] - 1) = *i;
        // reduce the number of elements we need to place in a a given area
        count[val] = count[val] - 1;
    }

}

/**
 * @brief Sorts 32-bit unsigned integers by calling 4 times the counting sort
 *         by sorting successvely 8-bit groups
 *
 * @note https://en.wikipedia.org/wiki/Radix_sort
 *
 * @param v vector to sort, modified by this function
 */
void RadixSort(std::vector<unsigned int> &v) {

    std::vector<unsigned> sorted(v.size(), 0);

    // loop through unsigned each bytes of unsigned ints
    for (size_t i = 0; i < sizeof(unsigned); ++i) {

        // create new Binary functor with the position
        //  of the current byte we are working on
        Binary<unsigned> binary(i);

        // sort the given vector using the Binary functor to sort
        //  using the Binary value.
        CountingSort(v.begin(), v.end(), sorted.begin(), binary);

        // save the current state of the sort
        v = sorted;
    }
}

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

// Ca fait le test 1 : des	séries	de	données	de	taille	différente de	même	distribution (meme vectors pour tous les tris)
// pa encore fini
void test1() {

    vector<int> vector_sizes {10, 100, 1000, 10000, 100000, 1000000};
    int randMin = 1;
    int randMax = 100;

    for (auto i = vector_sizes.begin(); i < vector_sizes.end(); ++i) {
        // Create a new vector contains : 1 ... currentSize)
        vector<int> v1(*i);
        for (int j = 0; j < *i; ++j) {
            int randNum = rand()%(randMax-randMin + 1) + randMin;
            v1.at(j) = randNum;
        }
        vector<int> v2 = v1;
        vector<int> v3 = v1;

        selectionSort(v1.begin(), v1.end());
        quickSort(v2.begin(), v2.end());
        // Add counting sort
    }

}


int main(int argc, const char * argv[]) {

    test1();

    return 0;
}