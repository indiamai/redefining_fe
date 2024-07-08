Examples
===========================

This page will describe the foundational element definitons and how they are written in this software.

DG0  
--------------------------
.. literalinclude:: ../../test/test_2d_examples_docs.py
   :language: python3
   :dedent:
   :start-after: [test_dg0 0]
   :end-before: [test_dg0 1]

DG1 on interval 
--------------------------

.. literalinclude:: ../../test/test_2d_examples_docs.py
   :language: python3
   :dedent:
   :start-after: [test_dg1_int 0]
   :end-before: [test_dg1_int 1]

DG1 on triangle
--------------------------
.. literalinclude:: ../../test/test_2d_examples_docs.py
   :language: python3
   :dedent:
   :start-after: [test_dg1_tri 0]
   :end-before: [test_dg1_tri 1]

CG elements
--------------------------

First, we set up the cells. As we are including our lower dim definitons into higher dim ones they need to resolve to common objects:

.. literalinclude:: ../../test/test_2d_examples_docs.py
   :language: python3
   :dedent:
   :start-after: [test_cg1 0]
   :end-before: [test_cg1 1]


Then, CG1 looks like this:

.. literalinclude:: ../../test/test_2d_examples_docs.py
   :language: python3
   :dedent:
   :start-after: [test_cg1 2]
   :end-before: [test_cg1 3]

And CG3 like this:

.. literalinclude:: ../../test/test_2d_examples_docs.py
   :language: python3
   :dedent:
   :start-after: [test_cg3 0]
   :end-before: [test_cg3 1]