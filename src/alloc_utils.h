#ifndef _ALLOC_UTILS
#define _ALLOC_UTILS

/*
   The alloctable is a chained list, where elemnts are added at the end
   and head is following the last element; next indicates the previous
   list element

   /head
   O<---O<---O<---O

   They are used to keep track of allocated arrays and free them.
 */

typedef struct alloctable {
    char name;			/* node identifier (for debugging)             */
    double *allocated;		/* Pointer to allocated array to keep track of */
    int keep_on_success;	/* If code executes cleanly, do not free when=1 */
    struct alloctable *next;	/* Pointer to previous node in list            */
} alloctable;

/* Add a node at the end of the list and move head to it */
void alloctable_add(alloctable ** head, double *allocated,
		    int keep_on_success, char name);

/* Free all allocated array in list but for keep_on_success flagged */
void alloctable_free_onsuccess(alloctable ** head);

/* Free the single node pointed by head and move it to next item */
void alloctable_free_last(alloctable ** head);

/* Free all allocated arrays in list */
void alloctable_free(alloctable ** head);

#endif
