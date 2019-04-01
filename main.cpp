/**
 * HEADER todo
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

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
 * @param output   beginning of the container where the result will be
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
        // if not, loop through to find it
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

    // loop through the count vector to set the range in which elements will be placed̉
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

void test1() {

    vector<unsigned> vectorSizes {10, 100, 1000, 10000, 100000, 1000000};
    const unsigned randMin = 1;
    const unsigned randMax = 100;
    const unsigned REPLICATION = 15;
    const double DIVISOR_NANO_TO_MILLIS = 1e+6;

    high_resolution_clock::time_point t1_selectionSort, t1_quickSort, t1_countingSort;
    high_resolution_clock::time_point t2_selectionSort, t2_quickSort, t2_countingSort;
    double selectionSortAverageTime = 0., quickSortAverageTime = 0., countingSortAverageTime = 0.;

    uniform_int_distribution<unsigned> alea (randMin, randMax);
    mt19937_64 gen(0);

    for (unsigned int &size : vectorSizes) {

        for (unsigned k = 0; k < REPLICATION; ++k) {

            // Create a new vector where it's size varies depending on the values in vectorSizes
            vector<unsigned> v1(size), v2(size), v3(size), w(size);

            // Generate random data in the vectors
            generate(v1.begin(), v1.end(), [&](){ return alea(gen); });
            v2 = v3 = v1;

            // Stop testing Selction Sort when we have vector sizes greater than 10'000
            //  it takes way too long after that
            if (size <= 10000) {
                t1_selectionSort = high_resolution_clock::now();
                selectionSort(v1.begin(), v1.end());
                t2_selectionSort = high_resolution_clock::now();
                selectionSortAverageTime += duration_cast<nanoseconds>(t2_selectionSort - t1_selectionSort).count();
            }


            // Calculate Quick sort execution time
            t1_quickSort = high_resolution_clock::now();
            quickSort(v2.begin(), v2.end());
            t2_quickSort = high_resolution_clock::now();
            quickSortAverageTime += duration_cast<nanoseconds>(t2_quickSort - t1_quickSort).count();


            // Calculate Counting sort execution time
            t1_countingSort = high_resolution_clock::now();
            CountingSort(v3.begin(), v3.end(), w.begin(), [&](unsigned value) {
                return value;
            }, alea.max());
            t2_countingSort = high_resolution_clock::now();
            countingSortAverageTime += duration_cast<nanoseconds>(t2_countingSort - t1_countingSort).count();

        }
        double finalTime;

        // Display the average time of execution for the Selection sort
        finalTime = selectionSortAverageTime / REPLICATION;
        if (size <= 10000) {
            cout << "For a vector of " << size << " selectionSort took "
                 << finalTime << " ns -> " << finalTime / DIVISOR_NANO_TO_MILLIS << " ms." << endl;
        }

        // Display the average time of execution for the Quick sort
        finalTime = quickSortAverageTime / REPLICATION;
        cout << "For a vector of " << size << " quickSort took "
             << finalTime << " ns -> " << finalTime / DIVISOR_NANO_TO_MILLIS << " ms." << endl;

        // Display the average time of execution for the Counting sort
        finalTime = countingSortAverageTime / REPLICATION;
        cout << "For a vector of " << size << " CountingSort took "
             << finalTime << " ns -> " << finalTime / DIVISOR_NANO_TO_MILLIS << " ms." << endl;

        cout << endl;
    }
}

void test2() {
    const unsigned vectorSize = 1000000;
    const unsigned minValue = 1;
    const vector<unsigned> maxValues {10, 100, 1000, 10000, 100000, 1000000};
    const unsigned REPLICATION = 15;
    const double DIVISOR_NANO_TO_MILLIS = 1e+6;

    high_resolution_clock::time_point t1_selectionSort, t1_quickSort, t1_countingSort;
    high_resolution_clock::time_point t2_selectionSort, t2_quickSort, t2_countingSort;
    double selectionSortAverageTime = 0., quickSortAverageTime = 0., countingSortAverageTime = 0.;

    vector<unsigned> v1(vectorSize), v2(vectorSize), v3(vectorSize);

    mt19937_64 gen(0);

    for (auto i = maxValues.begin(); i != maxValues.end(); ++i) {

        uniform_int_distribution<unsigned> alea (minValue, *i);

        for (unsigned k = 0; k < REPLICATION; ++k) {

            generate(v1.begin(), v1.end(), [&](){ return alea(gen); });

            v2 = v3 = v1;

            // quickSort time calcul
            t1_quickSort = high_resolution_clock::now();
            quickSort(v2.begin(), v2.end());
            t2_quickSort = high_resolution_clock::now();
            quickSortAverageTime += duration_cast<nanoseconds>(t2_quickSort - t1_quickSort).count();


            t1_countingSort = high_resolution_clock::now();
            RadixSort(v3);
            t2_countingSort = high_resolution_clock::now();
            countingSortAverageTime += duration_cast<nanoseconds>(t2_countingSort - t1_countingSort).count();

        }

        double finalTime;

        // QuickSort display average time
        finalTime = quickSortAverageTime / REPLICATION;
        cout << "For a vector of " << vectorSize << " [" << minValue << "," << (*i) <<"]" << " quickSort took "
             << finalTime << " ns -> " << finalTime / DIVISOR_NANO_TO_MILLIS << " ms." << endl;

        finalTime = countingSortAverageTime / REPLICATION;
        cout << "For a vector of " << vectorSize << " [" << minValue << "," << (*i) <<"]" << " RadixSort took "
             << finalTime << " ns -> " << finalTime / DIVISOR_NANO_TO_MILLIS << " ms." << endl;

        cout << endl;
    }
}

int main(int argc, const char * argv[]) {

    test2();

    return 0;
}