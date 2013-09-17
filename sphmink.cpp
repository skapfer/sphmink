#include "readpoly.h"
#include "Vector.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_specfunc.h>

#define MAX_L 12

static
double sign (double x)
{
    if (x<0.)
        return -1.;
    else
        return 1.;
}

class SphericalMinkowskis
{
    gsl_complex d_[91];
    double total_area;

    const gsl_complex *d (int l) const {
        return const_cast <SphericalMinkowskis *> (this)->d (l);
    }

    gsl_complex *d (int l) {
        switch (l)
        {
        case 0:
            return d_;
        case 1:
            return d_+1;
        case 2:
            return d_+3;
        case 3:
            return d_+6;
        case 4:
            return d_+10;
        case 5:
            return d_+15;
        case 6:
            return d_+21;
        case 7:
            return d_+28;
        case 8:
            return d_+36;
        case 9:
            return d_+45;
        case 10:
            return d_+55;
        case 11:
            return d_+66;
        case 12:
            return d_+78;
        case 13:
            return d_+91;
        default:
            std::cerr << "MAX_L broken\n";
            std::abort ();
        }
    }

public:
    SphericalMinkowskis () {
        std::memset (d_, 0, sizeof d_);
        total_area = 0.;
    }

    void add_facet (Vector f)
    {
        double area = norm (f);
        f /= area;
        double cos_th = f[2];
        double phi = atan2 (f[1], f[0]);
        total_area += area;
        gsl_complex *dp = d (0);
        for (int l = 0; l <= MAX_L; ++l)
        {
            double l_prefactor = sqrt (4*M_PI/(2*l+1));
            assert (d (l) == dp);
            for (int m = 0; m <= l; ++m, ++dp)
            {
                double leg = gsl_sf_legendre_sphPlm (l, m, cos_th);
                leg *= area;
                leg *= l_prefactor;
                gsl_complex ylm = gsl_complex_polar (leg, m*phi);
                *dp = gsl_complex_add (*dp, ylm);
            }
        }
    }

    double ql (int l) const
    {
        assert (l <= MAX_L);
        assert (l >= 0);
        const gsl_complex *pe = d (l+1);
        const gsl_complex *p = d (l);
        double r = gsl_complex_abs2 (*p);
        for (++p; p != pe; ++p)
            r += 2 * gsl_complex_abs2 (*p);
        return sqrt (r) / total_area;
    }

    double wl (int l) const
    {
        assert (l <= MAX_L);
        assert (l >= 0);
        gsl_complex v = {0., 0.};
        const gsl_complex *p = d (l);
        int ma, mb, mc;
        for (ma = -l; ma <= l; ++ma)
        for (mb = -l; mb <= l; ++mb)
        for (mc = -l; mc <= l; ++mc)
        {
            gsl_complex q = gsl_complex_polar (gsl_sf_coupling_3j (2*l, 2*l, 2*l, 2*ma, 2*mb, 2*mc), 0);

            // get the proper Clm coefficients from the arrays.
            // this is complicated a bit by the fact that we don't store the negative-m
            // coefficients, which are complex conjugates of the +m ones.
            gsl_complex a;
            if (ma < 0) {
                a = p[-ma];
                a = gsl_complex_conjugate (a);
                if ((-ma) % 2)
                    a = gsl_complex_mul_real (a, -1);
            }
            else
            {
                a = p[ma];
            }
            q = gsl_complex_mul (q, a);

            if (mb < 0) {
                a = p[-mb];
                a = gsl_complex_conjugate (a);
                if ((-mb) % 2)
                    a = gsl_complex_mul_real (a, -1);
            }
            else
            {
                a = p[mb];
            }
            q = gsl_complex_mul (q, a);

            if (mc < 0) {
                a = p[-mc];
                a = gsl_complex_conjugate (a);
                if ((-mc) % 2)
                    a = gsl_complex_mul_real (a, -1);
            }
            else
            {
                a = p[mc];
            }
            q = gsl_complex_mul (q, a);

            v = gsl_complex_add (v, q);
        }

        if (fabs (GSL_IMAG (v)) > .1e-4)
            std::cerr << "large spurious imaginary component (l = "
                      << l << GSL_REAL (v) << ", " << GSL_IMAG (v) << "i)\n";
        return pow (gsl_complex_abs (v), 1./3) / total_area * sign (GSL_REAL (v));
    }
};

class NormalDistribution : public PolyFileSink {
public:
    virtual
     vertex_id_t insert_vertex (double x, double y, double z,
                                const prop_list_t &)
    {
        long ret = vert_.size ();
        vert_.push_back (Vector (x, y, z));
        return ret;
    }

    virtual void insert_facet  (const vertex_id_list_t &vl,
                                const prop_list_t &props)
    {
        if (! (vl.size () >= 3))
            std::abort ();
        long label = props.label_from_alpha_value ();
        unsigned i;
        Vector area (0., 0., 0.);
        for (i = 0; i != vl.size () - 1; ++i)
            area += cross_product (getvert_ (vl[i]), getvert_ (vl[i+1]));
        area += cross_product (getvert_ (vl[i]), getvert_ (vl[0]));
        area /= 6;
        data[label].add_facet (area);
    }

    virtual void reexamine_vertex (vertex_id_t id, double *x) const
    {
        const Vector &y = getvert_ (id);
        *x++ = y[0]; *x++ = y[1]; *x = y[2];
    }

    void finalize ()
    {
        vert_.clear ();
    }

    std::map <long, SphericalMinkowskis> data;

private:
    const Vector &getvert_ (vertex_id_t id) const
    {
        return vert_.at (id);
    }
    std::vector <Vector> vert_;
};

int main (int, const char **argv)
{
    NormalDistribution nd;
    std::ifstream is (argv[1]);
    parse_poly_file (&nd, is);
    nd.finalize ();
    std::map <long, SphericalMinkowskis>::const_iterator it;
    for (it = nd.data.begin (); it != nd.data.end (); ++it)
    {
        std::cout << it->first << " ";
        for (int l = 0; l <= MAX_L; ++l)
            std::cout << it->second.ql (l) << " "
                      << it->second.wl (l) << " ";
        std::cout << "\n";
    }
    return 0;
}
