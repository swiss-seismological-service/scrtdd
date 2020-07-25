/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/


#ifndef __RTDD_APPLICATIONS_APP_H__
#define __RTDD_APPLICATIONS_APP_H__

#include <seiscomp3/client/streamapplication.h>

#include <iostream>
#include <vector>
#include <list>


class Application : public Seiscomp::Client::StreamApplication {
	public:
		DEFINE_SMARTPOINTER(Option);
		struct Option : public Seiscomp::Core::BaseObject {
			Option(const char *cfgname,
				   const char *cligroup = nullptr,
				   const char *cliparam = nullptr,
				   const char *clidesc = nullptr,
				   bool clidefault = false,
				   bool cliswitch = false)
			: cfgName(cfgname), cliGroup(cligroup),
			  cliParam(cliparam), cliDesc(clidesc),
			  cliDefault(clidefault), cliSwitch(cliswitch) {}

			virtual void bind(Seiscomp::System::CommandLine *cli) = 0;
			virtual bool get(Seiscomp::System::CommandLine *cli) = 0;
			virtual bool get(const Seiscomp::Client::Application *app) = 0;
			virtual void printStorage(std::ostream &os) = 0;

			const char *cfgName;
			const char *cliGroup;
			const char *cliParam;
			const char *cliDesc;
			bool cliDefault;
			bool cliSwitch;
		};

		typedef std::list<OptionPtr> Options;


	public:
		Application(int argc, char **argv);


	protected:
		void createCommandLineDescription();
		bool validateParameters();
		bool initConfiguration();

		void addOption(OptionPtr);

		void addOption(int *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(double *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(bool *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(std::string *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(std::vector<int> *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(std::vector<double> *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(std::vector<bool> *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		void addOption(std::vector<std::string> *var, const char *cfgname,
		               const char *cligroup = nullptr, const char *cliparam = nullptr,
		               const char *clidesc = nullptr, bool clidefault = false,
		               bool cliswitch = false);

		const Options &options() const;


	private:
		Options _options;
};


#endif
