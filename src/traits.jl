# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

divide(geotable::AbstractGeoTable) = values(geotable), domain(geotable)
attach(table, dom::Domain) = georef(table, dom)
