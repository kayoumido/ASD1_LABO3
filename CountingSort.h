//
//  CountingSort.h
//

#ifndef CountingSort_h
#define CountingSort_h

#include <vector>
#include <algorithm>

namespace asd1 {


    template<typename T>
    class Binary {
    public:
        Binary(size_t shift) : shift(shift) {}

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
}

#endif /* CountingSort_h */