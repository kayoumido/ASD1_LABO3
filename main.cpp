
#include <iostream>

// selectionSort
//
// Effectue le tri par sélection des éléments entre begin
// et end (non inclus). Doit appeler display() après chaque
// échange.
//
// A COMPLETER

template < typename RandomAccessIterator >
void selectionSort( RandomAccessIterator begin,
                    RandomAccessIterator end )
{

    for (RandomAccessIterator i = begin; i < end - 1; ++i) {
        RandomAccessIterator imin = i;

        for (RandomAccessIterator j = i + 1; j < end; ++j) {

            if (*j < *imin) {
                imin = j;
            }
        }

        swap(*i, *imin);
        display(begin, i, imin, end);
    }
}


int main(int argc, const char * argv[]) {
    return 0;

}