//This subroutine displaces the system to a new configuration
void mcmove(double x[], double y[], double z[], double diam[], double zz[], double *ener, int *nattemp, int *nacc, double del, double boxl) {
    *nattemp = *nattemp + 1;
    int no = (int)(ranf(&iseed) * np); // select random particle

    double xo = x[no];
    double yo = y[no];
    double zo = z[no];

    double r_cut = 0.2 * boxl;
    
    // Neighbor search
    int neighbors[np];
    int neighbor_count = 0;
    for (int j = 0; j < np; j++) {
        if (j == no) continue;
        double dx = x[j] - xo;
        double dy = y[j] - yo;
        double dz = z[j] - zo;
        double r2 = dx*dx + dy*dy + dz*dz;
        if (sqrt(r2) <= r_cut) {
            neighbors[neighbor_count++] = j;
        }
    }

    if (neighbor_count == 0) return; // no valid neighbor

    int nj = neighbors[rand() % neighbor_count]; // pick random neighbor

    // Compute original distance
    double dx_old = x[nj] - xo;
    double dy_old = y[nj] - yo;
    double dz_old = z[nj] - zo;
    double r12_old = sqrt(dx_old*dx_old + dy_old*dy_old + dz_old*dz_old);

    // Compute old energy only for moved particle
    double enero;
    energy(x, y, z, x[nj], y[nj], z[nj], diam, zz, &enero, nj);

    // Propose new configuration for nj
    double sigmaij = 0.5 * (diam[no] + diam[nj]);
    double r12_new = sigmaij + (r_cut - sigmaij) * ranf(&iseed);
    Vector3D dir = random_unit_vector(&iseed);

    double xj_old = x[nj];
    double yj_old = y[nj];
    double zj_old = z[nj];

    x[nj] = xo + r12_new * dir.x;
    y[nj] = yo + r12_new * dir.y;
    z[nj] = zo + r12_new * dir.z;

    // Apply periodic boundary conditions
    x[nj] -= boxl * round(x[nj] / boxl);
    y[nj] -= boxl * round(y[nj] / boxl);
    z[nj] -= boxl * round(z[nj] / boxl);

    // Reject immediately if overlap
    if (overlaps(nj, diam, x, y, z)) {
        x[nj] = xj_old;
        y[nj] = yj_old;
        z[nj] = zj_old;
        return;
    }

    // Compute new energy only for moved particle
    double enern;
    energy(x, y, z, x[nj], y[nj], z[nj], diam, zz, &enern, nj);

    // Metropolis criterion
    double deltaE = enern - enero;
    double acceptance = pow(r12_new / r12_old, 2) * exp(-beta * deltaE);
    
    if (ranf(&iseed) < acceptance) {
        *ener += 0.5 * deltaE; // update energy
        *nacc += 1;
    } else {
        // Revert move
        x[nj] = xj_old;
        y[nj] = yj_old;
        z[nj] = zj_old;
    }
}
