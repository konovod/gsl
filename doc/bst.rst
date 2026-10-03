.. index::
   single: binary search trees
   single: AVL trees
   single: red-black trees

*******************
Binary Search Trees
*******************

.. note::

   This chapter was generated automatically (with AI assistance) as a
   first draft and has not been reviewed by a human.  Check it against
   :file:`gsl_bst.h` and the sources under :file:`bst/` before relying on
   it.

The functions described in this chapter create and manipulate binary
search trees.  A binary search tree stores a set of items ordered by a
user-supplied comparison function, and supports insertion, lookup and
removal.  The module is derived from Ben Pfaff's GNU libavl and provides
AVL trees and red-black trees, both of which keep the height of the tree
logarithmic in the number of items.

The functions are declared in the header files :file:`gsl_bst.h`,
:file:`gsl_bst_types.h`, :file:`gsl_bst_avl.h` and :file:`gsl_bst_rb.h`.

Types
=====

.. type:: gsl_bst_type

   A structure containing the function pointers and other properties that
   define one kind of tree.  The library provides two instances,
   :data:`gsl_bst_avl` and :data:`gsl_bst_rb`; the routines below take one
   of these as their :data:`T` argument.

.. data:: const gsl_bst_type * gsl_bst_avl
          const gsl_bst_type * gsl_bst_rb

   The AVL tree and red-black tree types supplied by the library, for
   passing to :func:`gsl_bst_alloc`.

.. type:: gsl_bst_cmp_function

   The type of the comparison function::

      int compare (const void * a, const void * b, void * params);

   It must return a negative value if :data:`a` orders before :data:`b`,
   zero if the two are equal, and a positive value otherwise.

.. type:: gsl_bst_allocator

   A replacement memory allocator::

      typedef struct
      {
        void * (*alloc) (size_t size, void * params);
        void (*free) (void * block, void * params);
      } gsl_bst_allocator;

   A null pointer passed to :func:`gsl_bst_alloc` selects the standard
   :func:`malloc` and :func:`free` functions.

.. type:: gsl_bst_workspace

   A workspace holding a tree and its associated state.  It is created by
   :func:`gsl_bst_alloc` and released with :func:`gsl_bst_free`.

.. type:: gsl_bst_trav

   A traverser, used to walk through a tree in order.

Allocation
==========

.. function:: gsl_bst_workspace * gsl_bst_alloc (const gsl_bst_type * T, const gsl_bst_allocator * allocator, gsl_bst_cmp_function * compare, void * params)

   This function allocates and initializes a tree of type :data:`T`, which
   must be :data:`gsl_bst_avl` or :data:`gsl_bst_rb`.  The :data:`compare`
   function defines the ordering, and :data:`params` is passed unchanged
   to it and to the allocator.  On success the new workspace is returned;
   on failure a null pointer is returned and the error handler is invoked
   with :macro:`GSL_ENOMEM`.

.. function:: void gsl_bst_free (gsl_bst_workspace * w)

   This function frees every node in :data:`w` and the workspace itself.
   The items stored in the tree are not freed; if the tree owns them, the
   caller must free them first.

Modifying the tree
==================

.. function:: int gsl_bst_empty (gsl_bst_workspace * w)

   This function removes every item from :data:`w`, leaving it empty.

.. function:: void * gsl_bst_insert (void * item, gsl_bst_workspace * w)

   This function inserts :data:`item` into the tree.  If an equal item is
   already present, no insertion takes place and a pointer to the existing
   item is returned; otherwise :data:`item` is stored and returned.  A
   null pointer is returned if memory for the new node cannot be obtained.
   The library stores the pointer it is given and does not copy the item,
   so the item must stay valid while it is in the tree.

.. function:: void * gsl_bst_find (const void * item, const gsl_bst_workspace * w)

   This function returns a pointer to the item in the tree equal to
   :data:`item`, or a null pointer if there is none.

.. function:: void * gsl_bst_remove (const void * item, gsl_bst_workspace * w)

   This function removes the item equal to :data:`item` and returns a
   pointer to it.  If no matching item is present the tree is unchanged and
   a null pointer is returned.

.. function:: size_t gsl_bst_nodes (const gsl_bst_workspace * w)

   This function returns the number of items in the tree.

.. function:: size_t gsl_bst_node_size (const gsl_bst_workspace * w)

   This function returns the size in bytes of each node.

.. function:: const char * gsl_bst_name (const gsl_bst_workspace * w)

   This function returns a string naming the tree type, such as
   :code:`"AVL"` or :code:`"RB"`.

Traversal
=========

A traverser walks the items of a tree in their sorted order.  It must be
initialized by :func:`gsl_bst_trav_first`, :func:`gsl_bst_trav_last`,
:func:`gsl_bst_trav_find` or :func:`gsl_bst_trav_init` before use, and it
is invalidated by any modification of the tree.

.. function:: int gsl_bst_trav_init (gsl_bst_trav * trav, const gsl_bst_workspace * w)

   This function initializes :data:`trav` for the tree :data:`w`.  The
   traverser is not positioned on any item.

.. function:: void * gsl_bst_trav_first (gsl_bst_trav * trav, const gsl_bst_workspace * w)

   This function initializes :data:`trav` and returns the least-valued
   item in the tree, or a null pointer if the tree is empty.

.. function:: void * gsl_bst_trav_last (gsl_bst_trav * trav, const gsl_bst_workspace * w)

   This function initializes :data:`trav` and returns the greatest-valued
   item in the tree, or a null pointer if the tree is empty.

.. function:: void * gsl_bst_trav_find (const void * item, gsl_bst_trav * trav, const gsl_bst_workspace * w)

   This function searches for :data:`item`, positions :data:`trav` on it
   and returns it, or returns a null pointer if the item is not present.

.. function:: void * gsl_bst_trav_insert (void * item, gsl_bst_trav * trav, gsl_bst_workspace * w)

   This function behaves as :func:`gsl_bst_insert`, but also positions
   :data:`trav` on the inserted item (or on the existing duplicate).

.. function:: void * gsl_bst_trav_copy (gsl_bst_trav * dest, const gsl_bst_trav * src)

   This function copies the traverser :data:`src` into :data:`dest` and
   returns a pointer to the current item.

.. function:: void * gsl_bst_trav_next (gsl_bst_trav * trav)

   This function advances :data:`trav` to the next item in order and
   returns it, or returns a null pointer if there is none.

.. function:: void * gsl_bst_trav_prev (gsl_bst_trav * trav)

   This function moves :data:`trav` to the previous item in order and
   returns it, or returns a null pointer if there is none.

.. function:: void * gsl_bst_trav_cur (const gsl_bst_trav * trav)

   This function returns the current item without moving the traverser.

.. function:: void * gsl_bst_trav_replace (gsl_bst_trav * trav, void * new_item)

   This function replaces the current item with :data:`new_item` and
   returns the item it replaced.  The replacement must compare equal to the
   item it replaces, so that the ordering of the tree is unchanged.

References and Further Reading
==============================

The algorithms are those of the GNU libavl library:

* Ben Pfaff, *An Introduction to Binary Search Trees and Balanced Trees*,
  Free Software Foundation, 2004.
