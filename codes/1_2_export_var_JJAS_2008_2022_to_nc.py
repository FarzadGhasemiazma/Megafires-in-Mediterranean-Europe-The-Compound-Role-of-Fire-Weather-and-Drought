#!/usr/bin/env python

import xarray as xr
from pathlib import Path
import dask
import argparse

# Usar scheduler simple para no petar memoria
dask.config.set(scheduler="single-threaded")

MONTHS = [6, 7, 8, 9]  # JJAS
DATA_START = "2008-01-01"
DATA_END   = "2022-12-31"


def main():
    parser = argparse.ArgumentParser(
        description="Extrae VAR JJAS 2008–2022 de un Zarr y lo guarda en NetCDF."
    )
    parser.add_argument("--var", required=True, help="Nombre de la variable en el Zarr (p.ej. lai, tp, t2m)")
    parser.add_argument("--store", default="/diskonfire/MESOGEOS", help="Ruta al Zarr store")
    parser.add_argument("--outdir", default="/diskonfire/MESOGEOS", help="Directorio de salida para el NetCDF")
    args = parser.parse_args()

    var_name   = args.var
    store_path = args.store
    out_dir    = Path(args.outdir)
    out_nc     = out_dir / f"{var_name}_JJAS_2008_2022.nc"

    print(f"Abriendo Zarr: {store_path}")
    ds = xr.open_zarr(store_path, consolidated=False)

    if var_name not in ds:
        raise ValueError(f"'{var_name}' no está en el dataset. Variables: {list(ds.data_vars)}")

    da = ds[var_name]
    print(f"Variable seleccionada: {var_name} con dims {da.dims}")

    # 1) Seleccionar periodo 2008–2022
    da = da.sel(time=slice(DATA_START, DATA_END))

    # 2) Quedarse solo con meses JJAS (6,7,8,9)
    da = da.where(da.time.dt.month.isin(MONTHS), drop=True)

    print("Shape", var_name, "JJAS 2008–2022:", da.sizes)

    # 3) Rechunk conservador
    da = da.chunk({"time": 30, "y": 128, "x": 128})

    # 4) Encoding para NetCDF
    encoding = {
        var_name: {
            "zlib": True,
            "complevel": 3,
            "dtype": "float32",
            "chunksizes": (30, 128, 128),
        }
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    print("Escribiendo NetCDF en:", out_nc)

    da.to_dataset(name=var_name).to_netcdf(out_nc.as_posix(), encoding=encoding)

    print("Listo. NetCDF escrito en:", out_nc)


if __name__ == "__main__":
    main()
