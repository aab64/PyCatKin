.. _overview:
.. index:: Code overview


Code overview
*********************

Here is an overview of the energetic and kinetic equations and the modular structure of the code. 

Free energy components
-----------------------

The total free energy, :math:`G_{\textsf{total}}`, of a state is approximated by summation of electronic, translational, vibrational, and rotational energy contributions, i.e.,

.. math::
   :nowrap:
   
    \begin{equation}
        G_{\textsf{total}}
        =
        G_{\textsf{elec}}
        +
        G_{\textsf{tran}}
        +
        G_{\textsf{rota}}
        +
        G_{\textsf{vibr}}
        \textsf{.}
    \end{equation}

The electronic contribution, :math:`G_{\textsf{elec}}`, is taken from DFT. The translational contribution, :math:`G_{\textsf{tran}}`, is

.. math::
   :nowrap:
   
    \begin{equation}
        G_{\textsf{tran}}
        =
        -RT\ln\left(
        \frac{V}{N_{\textsf{A}}}
        {\left(\frac{2\pi Mk_{\textsf{B}}T}{h^2}\right)}^{\frac{3}{2}}
        \right)
        \textsf{,}
    \end{equation}

where :math:`V` is the standard molar volume and :math:`M` is the molecular mass. The rotational contribution, :math:`G_{\textsf{rota}}`, depends on the geometry of the molecule,

.. math::
   :nowrap:
   
    \begin{equation}
        G_{\textsf{rota}}
        =
        \left\lbrace
        \begin{array}{ll}
        -RT\ln\left(
        \frac{8\pi^2k_{\textsf{B}}TI}{\sigma h^2}
        \right)
        &
        \textsf{linear molecule}
        \\
        -RT\ln\left(
        \frac{\sqrt{\pi}}{\sigma}
        {\left(\frac{8\pi^2k_{\textsf{B}}T}{h^2}\right)}^{\frac{3}{2}}
        \sqrt{I_1I_2I_3}
        \right)
        & \textsf{otherwise.}
        \end{array}
        \right.
    \end{equation}

Here, :math:`I` and :math:`I_1`, :math:`I_2`, :math:`I_3` are the inertia and three principle components of the inertia, respectively and :math:`\sigma` is the symmetry number of the molecule. The vibrational contribution is computed using the vibrational frequencies, :math:`v_i`, accounting for all degrees of freedom not attributed to overall translational or rotational movement,

.. math::
   :nowrap:
   
    \begin{equation}
        G_{\textsf{vibr}}
        =
        R\sum_{i=1}^{N_{\textsf{dof}}}
        \left[
        \frac{hv_i}{2k_{\textsf{B}}}
        +
        T\ln\left(
        1-e^{-\frac{hv_i}{k_{\textsf{B}}T}}
        \right)
        \right]
        \textsf{.}
    \end{equation}

Molecules interacting with the surface can be assumed to experience no/hindered/full translational or rotational contributions. Molecules in the gas-phase have three translational and two/three rotational contributions for linear/nonlinear molecules respectively. The number of atoms comprising a state includes all mobile surface atoms, :math:`N_{\textsf{s}}`, atoms in the molecules interacting with the surface, :math:`N_{\textsf{m}}`, and atoms in the molecules moving freely as gases, :math:`N_{\textsf{g}}`. Thus, the number of degrees of freedom is:

.. math::
   :nowrap:
   
    \begin{equation}
        N_{\textsf{dof}}
        =
        \left\lbrace
        \begin{array}{ll}
            3N_{\textsf{s}} + 3N_{\textsf{m}} & \textsf{Surface-bound state}\\
            3N_{\textsf{s}} + 3N_{\textsf{g}}  - 5/6 & \textsf{Surface and linear/nonlinear gas molecule}\\
            3N_{\textsf{g}} - 5/6 & \textsf{Linear/nonlinear gas molecule.}
        \end{array}
        \right.
    \end{equation}
    
The hindered model assumes that the surface-bound state has some fraction, :math:`\alpha`, of the translational and rotational energy of the corresponding free gas molecule:

.. math::
   :nowrap:
   
    \begin{equation}
        G_{\textsf{total}}
        =
        G_{\textsf{elec}}
        +
        \alpha\left(G_{\textsf{tran}}
        +
        G_{\textsf{rota}}\right)
        +
        G_{\textsf{vibr}}
        \textsf{.}
    \end{equation}

Reaction barriers and energy
------------------------------

The forward and reverse reaction barriers are computed as the difference in Gibbs free energy between the transition state (TS) and initial state (IS) for the forward reaction, 

.. math::
   :nowrap:
   
    \begin{equation}
        \Delta G_{\textsf{a}}^{\textsf{f}}
        =
        G_{\textsf{total,TS}} - G_{\textsf{total,IS}}
        \text{,}
    \end{equation}

and the transition state and the final state (FS) for the reverse reaction,

.. math::
   :nowrap:
   
    \begin{equation}
        \Delta G_{\textsf{a}}^{\textsf{r}}
        =
        G_{\textsf{total,TS}} - G_{\textsf{total,FS}}
        \text{,}
    \end{equation}

respectively. The reaction energies are computed as the difference in Gibbs free energy between each pair of initial and final states, 

.. math::
   :nowrap:
   
    \begin{equation}
        \Delta G_{\textsf{r}}
        =
        G_{\textsf{total,FS}} - G_{\textsf{total,IS}}
        \text{.}
    \end{equation}

Energy span model
---------------------

The energy span  (ES) model is a way of describing the theoretical efficiency of the catalytic cycles using state energies derived from the first-principles calculations (Kozuch and Shaik, 2011). The TOF is given by the summation of the pairwise energy differences,

.. math::
   :nowrap:
   
    \begin{equation}
        \textsf{TOF}
        =
        \frac{k_{\textsf{B}}T}{h}
        \frac{e^{-\Delta G_{\textsf{r}}/RT} - 1}
        {\sum_{ij} e^{\left(G_{T_i} - G_{I_j} - \delta G_{ij}\right)/RT}}
        \textsf{,}
    \end{equation}

where the summation is taken over all transition states (:math:`T`) and reaction intermediates (:math:`I`) in the catalytic cycle. Indexing states :math:`T` by :math:`i` and :math:`I` by :math:`j`, the term :math:`\delta G_{ij}` is defined as follows:

.. math::
   :nowrap:
    
    \begin{equation}
        \delta G_{ij}
        =
        \left\lbrace
        \begin{array}{cc}
            \Delta G_{\textsf{r}} & i\geq j\\
            0 & i< j\textsf{.}
        \end{array}
        \right.
    \end{equation}

The pairwise energy differences can be used to characterize the degree of TOF control of transition states,

.. math::
   :nowrap:
    
    \begin{equation}
        X_{\textsf{TOF,T}_i}
        =
        \frac{\sum_j e^{\left(G_{T_i} - G_{I_j} - \delta G_{ij}\right)/RT}}
        {\sum_{ij} e^{\left(G_{T_i} - G_{I_j} - \delta G_{ij}\right)/RT}}
        \textsf{,}
    \end{equation}

and intermediates,

.. math::
   :nowrap:
    
    \begin{equation}
        X_{\textsf{TOF,I}_j}
        =
        \frac{\sum_i e^{\left(G_{T_i} - G_{I_j} - \delta G_{ij}\right)/RT}}
        {\sum_{ij} e^{\left(G_{T_i} - G_{I_j} - \delta G_{ij}\right)/RT}}
        \textsf{.}
    \end{equation}

In this way, we can determine which states are the most TOF controlling, the so-called turnover-determining transition state (TDTS) and turnover-determining intermediate (TDI). These are the states with the largest energy difference across successive catalytic cycles.

Microkinetic model
---------------------

A mean-field microkinetic model describes elementary steps occuring on a catalyst surface, including adsorption/desorption and reactions, assuming that the surface coverage is homogeneous. Reaction steps are formulated using the law of mass action. For each elementary step, :math:`i`, the reaction rate is given by,

.. math::
   :nowrap:
   
    \begin{equation}
        r_i
        =
        k_i^{\textsf{f}}\prod_j\theta_{ij}\prod_j p_{ij}-
        k_i^{\textsf{r}}\prod_l\theta_{il}\prod_l p_{il}
    \end{equation}

and for each species, :math:`j`, the differential equation governing its rate of change is:

.. math::
   :nowrap:
   
    \begin{equation}
        \frac{\partial\theta_j}{\partial t}
        =
        \sum_i \nu_{ij}r_i
    \end{equation}

where :math:`\nu_{ij}` is the stoichiometry of species :math:`j` in reaction :math:`i`. In order to conserve the number of surface sites, with :math:`\theta^{\textsf{tot}}` the normalised number of sites,

.. math::
   :nowrap:
   
    \begin{equation}
        \sum_j\theta_j
        =
        \theta^{\textsf{tot}}
        \textsf{,}
    \end{equation}

if there is a single site type, or 

.. math::
   :nowrap:
   
    \begin{equation}
        \sum_{j\in S_k}\theta_j
        =
        \theta^{\textsf{tot}}_k
        \textsf{,}
    \end{equation}

if there are multiple site types :math:`S_k`. This is achieved naturally by tracking all clean surface site types in the same manner as adsorbates. Then, the Jacobian for this system of differential equations is:

.. math::
   :nowrap:
   
    \begin{equation}
        \frac{\partial}{\partial{\theta}_n}
        \frac{\partial{\theta}_j}{\partial t}
        =
        \sum_i \nu_{ij}\frac{\partial r_i}{\partial{\theta}_n}
        =
        \sum_i \nu_{ij}
        \left[
        k_i^{\text{f}}\prod_m p_{im}\sum_q
        \delta_{qn}
        \prod_{m\ne q}\theta_{im}-
        k_i^{\text{r}}\prod_l p_{il}\sum_s
        \delta_{sn}
        \prod_{l\ne s}\theta_{il}
        \right]
    \end{equation}

Reaction rate constants are typically taken to have an Arrhenius form, with the pre-factor determined from transition state theory,

.. math::
   :nowrap:
    
    \begin{equation}
        k
        =
        \frac{k_{\textsf{B}}T}{h}
        \exp\left(
        {-\frac{\Delta G_{\textsf{a}}}{RT}}
        \right)
        \textsf{,}
    \end{equation}

where :math:`k_{\textsf{B}}` is the Boltzmann constant, :math:`h` is the Planck constant, :math:`T` is the temperature, :math:`R` is the gas constant and :math:`\Delta G_{\textsf{a}}` is the activation free energy. The pre-factors for adsorption rate constants are instead determined using collision theory,

.. math::
   :nowrap:
    
    \begin{equation}
        k_{\textsf{ads}}
        =
        \frac{A}{\sqrt{2\pi M k_{\textsf{B}}T}}
        \textsf{,}
    \end{equation}

where :math:`A` is the area of a site and :math:`M` is the mass of the molecule. The reverse reaction rate constants can be computed from the corresponding equilibrium constants,

.. math::
   :nowrap:
    
    \begin{equation}
        K_{\textsf{eq}}
        =
        \exp\left(
        -\frac{\Delta G_{\textsf{r}}}{RT}
        \right)
        \textsf{.}
    \end{equation}

using the reaction free energy, :math:`\Delta G_{\textsf{r}}` to ensure thermodynamic consistency. 

The degree of rate control (DRC) is computed from the overall reaction rate, 

.. math::
   :nowrap:
   
    \begin{equation}
        \chi_i 
        = 
        {\left(
        \frac{\partial \ln r}
        {\partial \ln k_i}
        \right)}
        _{k_j\ne k_i, K_i}
        =
        \frac{k_i}{r}
        {\left(\frac{\partial r}{\partial k_i}\right)}_{k_j\ne k_i, K_i}
        \textsf{,}
    \end{equation}

and can be estimated numerically using finite differences to perturb the rate constants and describe the derivative. 

Reactor models
---------------------

Here, we consider a continuously stirred tank reactor (CSTR) model in which the mixture is assumed to be spatially homogeneous. A CSTR is parameterized by its residence time, :math:`\tau` (i.e., the reactor volume, :math:`V`, divided by the flow rate, :math:`Q`) and, for heterogeneous reactions occurring on a solid catalyst surface, the total number of catalyst sites (i.e., the site density, :math:`\rho_{\text{cat}}`, times the catalyst area, :math:`A_{\text{cat}}`):

.. math::
   :nowrap:
   
    \begin{equation}
        \frac{dp_{i}}{dt}
        =
        \tau^{-1}\left(p_{i}^{\textsf{in}} - p_{i}\right)
        +
        \frac{k_{\textsf{B}}T}{V} N_{\textsf{sites}} S_i
        \textsf{,}
    \end{equation}

where :math:`S_i` represents the sink/source term for mass transport to/from the surface from the gas due to adsorption/desorption respectively and 

.. math::
   :nowrap:
   
    \begin{equation}
        \tau = \frac{V}{Q}
        \textsf{,}
    \end{equation}

and

.. math::
   :nowrap:
   
    \begin{equation}
        N_{\textsf{sites}} = \rho_{\text{cat}}A_{\text{cat}}
        \textsf{.}
    \end{equation}

Structure of modules
---------------------
:program:`PyCatKin` is written using `object-oriented <https://docs.python.org/3/tutorial/classes.html>`_ programming.
The central subpackage modules are defined as shown in the figure below:

.. only:: latex
   
   .. image:: source/code_layout/code_layout.pdf
      :width: 750

.. only:: html
   
   .. image:: source/code_layout/code_layout.svg
      :width: 750

The modules have the following functions:
    - **State**: A microscopic state, which can be either a surface, adsorbant(s) or gas molecule(s). Information stored can include mass and inertia, energetic terms, structure.
    - **ScalingState**: A microscopic state, which can be either a surface, adsorbant(s) or gas molecule(s). Superficially the same as an instance of State; however, energetic terms are defined by scaling relation.
    - **Energy**: An ordered collection of states defining a reaction energy landscape (which can be drawn). Energy span model calculations can be performed using its member functions.
    - **Reaction**: An elementary surface reaction. Defined by lists of reactant, product and transition state states, the site area on which the reaction occurs and the scaling (for example due to a non-unity sticking coefficient). Used to compute reaction energies and barriers, and thus, reaction rate constants. 
    - **Reactor**: The system being studied. Can be either a CSTR or an infinite dilution reactor. Provides information about the boundary conditions (mass transport).
    - **System**: Defined by a set of states, a set of reactions, and a reactor, with parameters for initial and boundary conditions, and solver tolerances. Used to solve for transient or steady-state profiles of relevant species, compute the DRC, and to save results.
    - **Uncertainty**: Defined by a system, the mean and variance of Gaussian distribution, and the number of samples. Used to estimate the propagation of uncertainty in the energy landscape to the kinetics with a correlated error model. 

The reaction rate constants, input parser and some preset simulations are defined in the functions subpackage and the physical constants are defined in the constants subpackage. 
