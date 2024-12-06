/* Contains utilities to work with allocated arrays */
#include "alloc_utils.h"
#include <stdlib.h>

/*
   - head is the address of the pointer pointing to last element in list;
   - next_node needs to be malloc'd to make it persistent
   - next_node's next points to current head; then head is moved
 */
void alloctable_add(alloctable ** head, double *allocated,
		    int keep_on_success, char name)
{
    alloctable *next_node = malloc(1 * sizeof(alloctable));
    next_node->allocated = allocated;
    next_node->keep_on_success = keep_on_success;
    next_node->next = *head;
    next_node->name = name;
    *head = next_node;
}

/* Free allocated array and alloctable pointer */
void alloctable_free_last(alloctable ** head)
{
    alloctable *curr_node = *head;
    if (curr_node->allocated) {
//	printf("DEALLOCATING: %c\n", curr_node->name);
	free(curr_node->allocated);
    }
    *head = curr_node->next;
    free(curr_node);
}

void alloctable_free(alloctable ** head)
{
    while (*head)
	alloctable_free_last(head);
}



void alloctable_free_onsuccess(alloctable ** head)
{
    while (*head) {
	alloctable *curr_node = *head;
	if (curr_node->allocated && !curr_node->keep_on_success) {
//	    printf("DEALLOCATING: %c\n", curr_node->name);
	    free(curr_node->allocated);
	}
	*head = curr_node->next;
	free(curr_node);
    }
}
