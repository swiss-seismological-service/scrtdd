.. _database-label:

Database connection
===================

When SeisComP modules need to access the database for reading or writing data (events, picks, magnitudes, etc.) they use the connection string configured in either ``global.cfg`` (which is inherited by every module) or in ``scmaster.cfg``, in which case is scmaster module that passes the database connection string to every module when they connect to the messaging system (usually at module startup).

However, when running rtDD from the command line, it doesn't connect to the messaging system and if the database connection is specified via ``scmaster.cfg``, the information never reaches rtDD. In this case the database connection must be passed as a command line option::

    scrtdd [some options] -d  mysql://user:password@host/seiscompDbName

or in case of a Postgresql database::

    scrtdd [some options] --plugins dbpostgresql -d postgresql://user:password@host/seiscompDbName

It is worth noting that this feature allows rtDD to connect to remote seiscomp databases too and relocate events stored in other machines without interfering with the real-time processing happening there.

