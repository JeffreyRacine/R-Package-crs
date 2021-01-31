/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search -                */
/*                                                                                 */
/*  NOMAD - version 3.9.1 has been created by                                      */
/*                 Charles Audet               - Ecole Polytechnique de Montreal   */
/*                 Sebastien Le Digabel        - Ecole Polytechnique de Montreal   */
/*                 Viviane Rochon Montplaisir - Ecole Polytechnique de Montreal   */
/*                 Christophe Tribes           - Ecole Polytechnique de Montreal   */
/*                                                                                 */
/*  The copyright of NOMAD - version 3.9.1 is owned by                             */
/*                 Sebastien Le Digabel        - Ecole Polytechnique de Montreal   */
/*                 Viviane Rochon Montplaisir - Ecole Polytechnique de Montreal   */
/*                 Christophe Tribes           - Ecole Polytechnique de Montreal   */
/*                                                                                 */
/*  NOMAD v3 has been funded by AFOSR and Exxon Mobil.                             */
/*                                                                                 */
/*  NOMAD v3 is a new version of NOMAD v1 and v2. NOMAD v1 and v2 were created     */
/*  and developed by Mark Abramson, Charles Audet, Gilles Couture, and John E.     */
/*  Dennis Jr., and were funded by AFOSR and Exxon Mobil.                          */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Ecole Polytechnique de Montreal - GERAD                                      */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/

/**
 \file   XMesh.hpp
 \brief  Class for the xmesh (extensible mesh) (headers)
 \author Christophe Tribes
 \date   2014-07
 \see    XMesh.cpp
 */
#ifndef __XMESH__
#define __XMESH__

#include "OrthogonalMesh.hpp"

namespace NOMAD {
    
    /// Class for the MADS Xmesh.
    /**
     - The Xmesh in NOMAD is defined with the basic directions and the
     mesh size parameter delta^k for each coordinate.
     */
    
    class XMesh : public NOMAD::OrthogonalMesh {
        
        /*--------------------------------------------------------------*/
    private:
        
        NOMAD::Point    _r; // Mesh index per coordinate.
        NOMAD::Point    _r_min;
        NOMAD::Point    _r_max;
        
        /*--------------------------------------------------------------*/
        
        /// Private affectation operator.
        /**
         \param m The right-hand side object -- \b IN.
         */
        const XMesh & operator = ( const XMesh & m );
        
        /// Check the minimal poll size criterion.
        bool check_min_poll_size_criterion ( ) const;
        
        /// Check the minimal mesh size criterion.
        bool check_min_mesh_size_criterion ( ) const;
        
        
        /// Access to the mesh size parameter delta^k.
        /**
         \return delta    The mesh size parameter delta^k -- \b OUT.
         \param  i        The index of the mesh size
         */
        inline NOMAD::Double get_delta ( int i ) const ;
        
        
        /// Access to the poll size parameter Delta^k.
        /**
         \return Delta    The poll size parameter Delta^k -- \b OUT.
         \param  i        The index of the poll size
         */
        inline NOMAD::Double get_Delta ( int i ) const ;
        
        void init ( );
        
    public:
        
        /// Constructor.
        /**
         \param anisotropic_mesh        Use anisotropic mesh or not                     -- \b IN.
         \param Delta_0                 Initial poll size Delta^0                       -- \b IN.
         \param Delta_min               Minimal poll size Delta^min (may be undefined)  -- \b IN.
         \param delta_min               Minimal mesh size delta^min (may be undefined)  -- \b IN.
         \param poll_update_basis       Poll update basis (b); default=2                -- \b IN.
         \param poll_coarsening_step    Poll coarsening step (w+); default=1            -- \b IN.
         \param poll_refining_step      Poll refining step (w-); default=-1             -- \b IN.
         \param fixed_variables         Fixed variables                                 -- \b IN.
         \param limit_min_mesh_index    Limit mesh index (<0)                           -- \b IN.
         */
        XMesh (bool                 anisotropic_mesh,
               const NOMAD::Point & Delta_0   ,
               const NOMAD::Point & Delta_min ,
               const NOMAD::Point & delta_min  ,
               const NOMAD::Point & fixed_variables,
               NOMAD::Double        poll_update_basis=2.0,
               int                  poll_coarsening_step=1,
               int                  poll_refining_step=-1 ,
               int                  limit_min_mesh_index=NOMAD::XL_LIMITS)
        : NOMAD::OrthogonalMesh ( anisotropic_mesh,
                                 NOMAD::Double(),
                                 Delta_0,
                                 Delta_min,
                                 delta_min,
                                 fixed_variables,
                                 NOMAD::Point(),
                                 poll_update_basis,
                                 poll_coarsening_step,
                                 poll_refining_step ,
                                 limit_min_mesh_index ) { init();}
        
        
        
        
        /// Copy constructor.
        /**
         \param m The copied object -- \b IN.
         */
        XMesh ( const XMesh & m )
        : OrthogonalMesh ( m ) ,
        _r( m._r ),
        _r_min      ( m._r_min ),
        _r_max      ( m._r_max ) { init(); }
        
        /// Destructor.
        ~XMesh ( void )
        {
            _delta_0.clear();
            _Delta_0.clear();
            _delta_min.clear();
            _Delta_min.clear();
        }
        
        
        /// Access to the mesh size parameter delta^k for earch coordinate.
        /**
         - It is a NOMAD::Point of size \c nc the number of free variables.
         \param delta    The mesh size parameter delta^k -- \b OUT.
         \return A boolean equal to \c true if all values are
         strictly inferior than the associated minimal
         mesh size delta_min (stopping criterion MIN_MESH_SIZE).
         */
        bool get_delta ( NOMAD::Point & delta ) const ;
        
        
        /// Access to the largest mesh size for earch coordinate.
        /**
         \return delta_max
         */
        NOMAD::Point get_delta_max ( void ) const { return _delta_0;}
        
        
        /// Access to the poll size Delta^k for each coordinate.
        /**
         - It is a NOMAD::Point of size \c nc the number of free variables.
         \param Delta    The poll size parameter Delta^k -- \b OUT.
         \return A boolean equal to \c true if all values are
         strictly inferior than the associated minimal
         poll size Delta_min
         (stopping criterion MIN_POLL_SIZE).
         */
        bool get_Delta ( NOMAD::Point & Delta) const ;
        
        
        /// Update the XMesh (poll and mesh sizes).
        /**
         \param success    Type of success of the iteration     -- \b IN.
         \param dir        Direction of the iteration           -- \b IN.
         */
        void update ( NOMAD::success_type success, const NOMAD::Direction * dir=NULL);
        
        /// Reset the mesh to its original size (mesh indices).
        void reset ( void ) { init() ;}
        
        
        /// Display.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display ( const NOMAD::Display & out ) const;
        
        /// Test if r < r_min so far for all coordinates.
        /**
         \return True if mesh is the finest so far, False otherwise.
         */
        bool is_finest ( void ) const;
        
        
        /// Scale and project the ith component of a vector on the mesh
        /**
         \param i        The vector component number        -- \b IN.
         \param l        The vector component value         -- \b IN.
         \param round_up The flag for round-up              -- \b IN.
         \return         The ith component of a vector after mesh scaling and projection
         */
        NOMAD::Double scale_and_project(int i, const NOMAD::Double & l , bool round_up ) const ;
        
        
        /// Check the stopping conditions on the minimal poll and mesh sizes.
        /**
         \param stop           Stop flag                  -- \b IN/OUT.
         \param stop_reason    Stop reason                -- \b OUT.
         */
        void check_min_mesh_sizes (bool             & stop           ,
                                   NOMAD::stop_type & stop_reason      ) const;
        
        /// Access to the mesh indices per coordinate.
        /**
         \return A point with the mesh index for each coordinate.
         */
        const NOMAD::Point get_mesh_indices ( void  ) const { return _r; }
        
        /// Access to the min mesh indices reached so far.
        /**
         \return A point with the mesh index for each coordinate.
         */
        const NOMAD::Point get_min_mesh_indices ( void  ) const { return _r_min; }
        
        
        /// Access to the max mesh indices reached so far.
        /**
         \return A point with the mesh index for each coordinate.
         */
        const NOMAD::Point get_max_mesh_indices ( void  ) const { return _r_max; }
        
        
        /// Manually set the mesh index.
        /**
         \param r  The mesh index provided as a point -- \b IN.
         */
        void set_mesh_indices ( const NOMAD::Point & r );
        
        
        
        /// Manually set the limit mesh index (termination criterion).
        /**
         \param l  The limit mesh index for all coordinates -- \b IN.
         */
        void set_limit_mesh_index ( int l );
        
        
        
        /// Access to the mesh ratios after a success
        /**
         \return A point with the ratio for each coordinate
         */
        NOMAD::Point get_mesh_ratio_if_success ( void ) const;
        
        
        
    };
    
    /// Display a NOMAD::XMesh object.
    /**
     \param out The NOMAD::Display object -- \b IN.
     \param m   The NOMAD::XMesh object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
                                               const NOMAD::XMesh    & m     )
    {
        m.display ( out );
        return out;
    }
    
}

#endif
