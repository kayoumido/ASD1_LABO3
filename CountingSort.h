//
//  CountingSort.h
//

#ifndef CountingSort_h
#define CountingSort_h

#include <vector>
#include <algorithm>
#include <type_traits>

namespace asd1 {

    /**
     Tri comptage générique

     https://en.wikipedia.org/wiki/Counting_sort

     @param first  [first,last] est la plage d'éléments d'un tableau à trier.
     @param last   [first,last] est la plage d'éléments d'un tableau à trier.
     @param output début du tableau où écrire le résultat. Doit être différent
                   du tableau d'entrée
     @param key fonction prenant un élément en entrée et retourant sa position
                dans le tableau de comptage
     @param max_key valeur maximale pouvant être retournée par key(...). Si -1,
                    la fonction calcule cette valeur en parcourant le tableau une
                    fois de plus
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

        for (auto i = first; i != last; ++i) {
            auto val = key(*i);
            count[val] = count[val] + 1;
        }

        for (size_t i = 1; i < count.size(); ++i) {
            count[i] = count[i] + count[i - 1];
        }

        for (auto i = last; i != first; ) {
            --i;
            auto val = key(*i);

            *(output + count[val] - 1) = *i;
            count[val] = count[val] - 1;
        }

    }

    /**
     Tri par base d'entiers non signés sur 32 bits mis en oeuvre en appelant
     4 fois le tri comptage en triant successivement par groupe de 8 bits.

     https://en.wikipedia.org/wiki/Radix_sort

     @param v vecteur à trier, modifié par cette fonction
     */
    void RadixSort(std::vector<unsigned int> &v) {
// a compléter
    }
}

#endif /* CountingSort_h */