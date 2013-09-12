My First Projection
===================

Step 1.
-------

Solving the Poission equation of :math:`u^*` 

Step 2.
-------

Project the inter velocity into divergence space

.. code-block:: matlab

   function Ku = getMatrix(n,h)

	if isscalar(h)
    	h = h.*ones(n-1,1);
	end

	Dp = spdiags(h(2:end),0,n-1,n-1)\spdiags([-1,0;ones(n-3,1)*[-1,1];0,1],[0,1],n-1,n-1);
	Dn = spdiags(h(1:end-1),0,n-1,n-1)\spdiags([-1,1;ones(n-3,1)*[-1,1];0,1],[-1,0],n-1,n-1);

	Ku = spdiags(h(2:end)+h(1:end-1),0,n-1,n-1)\ (2*(Dp-Dn));
	end

Step 3.
-------

Update pressure