#!/usr/bin/env seiscomp-python
# -*- coding: utf-8 -*-
############################################################################
# Copyright (C) GFZ Potsdam, ETH Zuerich                                   #
#                                                                          #
# All rights reserved.                                                     #
#                                                                          #
# GNU Affero General Public License Usage                                  #
# This file may be used under the terms of the GNU Affero                  #
# Public License version 3.0 as published by the Free Software Foundation  #
# and appearing in the file LICENSE included in the packaging of this      #
# file. Please review the following information to ensure the GNU Affero   #
# Public License version 3.0 requirements will be met:                     #
# https://www.gnu.org/licenses/agpl-3.0.html.                              #
############################################################################

import sys
import seiscomp.core as sc_core
import seiscomp.client as sc_client
import seiscomp.datamodel as sc_datamodel
import seiscomp.logging as sc_logging
from collections import namedtuple


def _parseTime(timestring):
    t = sc_core.Time()
    if t.fromString(timestring, "%F %T"):
        return t
    if t.fromString(timestring, "%FT%T"):
        return t
    if t.fromString(timestring, "%FT%TZ"):
        return t
    return None


class EventList(sc_client.Application):

    def __init__(self, argc, argv):
        sc_client.Application.__init__(self, argc, argv)

        self.setMessagingEnabled(False)
        self.setDatabaseEnabled(True, False)
        self.setDaemonEnabled(False)

        self._startTime = None
        self._endTime = None
        self._modifiedAfterTime = None

    def createCommandLineDescription(self):
        self.commandline().addGroup("Events")
        self.commandline().addStringOption("Events", "begin",
                                           "specify the lower bound of the time interval e.g. 2000-01-01 00:00:00")
        self.commandline().addStringOption(
            "Events", "end", "specify the upper bound of the time interval e.g. 2000-01-01 00:00:00")
        self.commandline().addStringOption("Events", "modified-after",
                                           "select events modified after the specified time e.g. 2000-01-01 00:00:00")
        self.commandline().addStringOption("Events", "ev-type",
                                           "include only events whose type is one of the values provided (comma separated list) "
                                           'e.g. --ev-type "earthquake,induced earthquake. An empty string matches events with '
                                           'no type e.g --ev-type "" or --ev-type "earthquake,"')
        self.commandline().addOption(
            "Events", "simple", "Simple output. Print only origin ids")
        self.commandline().addGroup("Origins")
        self.commandline().addStringOption(
            "Origins", "org-type", "preferred, last, first, maxPhases, minPhases, maxRMS, minRMS (default is preferred)")
        self.commandline().addOption(
            "Origins", "manual-only", "Include only manual origins")
        self.commandline().addOption(
            "Origins", "auto-only", "Inlude only automatic origins")
        self.commandline().addStringOption("Origins", "inc-author",
                                           "include only origins whose author is one of the values provided (comma separated list)")
        self.commandline().addStringOption("Origins", "excl-author",
                                           "exclude origins whose author is one of the values provided (comma separated list)")
        self.commandline().addStringOption("Origins", "inc-method",
                                           "include only origins whose methodID is one of the values provided (comma separated list)")
        self.commandline().addStringOption("Origins", "excl-method",
                                           "exclue origins whose methodID is one of the values provided (comma separated list)")
        self.commandline().addStringOption("Origins", "inc-agency",
                                           "include only origins whose agencyID is one of the values provided (comma separated list)")
        self.commandline().addStringOption("Origins", "excl-agency",
                                           "exclude origins whose agencyID is one of the values provided (comma separated list)")
        self.commandline().addStringOption("Origins", "area",
                                           "Include only origins in the rectangular area provided: MinLat,MinLon,MaxLat,MaxLon")
        return True

    def init(self):
        if not sc_client.Application.init(self):
            return False

        try:
            start = self.commandline().optionString("begin")
        except BaseException:
            start = "1900-01-01T00:00:00Z"
        self._startTime = _parseTime(start)
        if self._startTime is None:
            sc_logging.error("Wrong 'begin' format '%s'" % start)
            return False
        sc_logging.debug("Setting start to %s" %
                         self._startTime.toString("%FT%TZ"))

        try:
            end = self.commandline().optionString("end")
        except BaseException:
            end = "2500-01-01T00:00:00Z"
        self._endTime = _parseTime(end)
        if self._endTime is None:
            sc_logging.error("Wrong 'end' format '%s'" % end)
            return False
        sc_logging.debug("Setting end to %s" %
                         self._endTime.toString("%FT%TZ"))

        try:
            modifiedAfter = self.commandline().optionString("modified-after")
            self._modifiedAfterTime = _parseTime(modifiedAfter)
            if self._modifiedAfterTime is None:
                sc_logging.error(
                    "Wrong 'modified-after' format '%s'" % modifiedAfter)
                return False
            sc_logging.debug(
                "Setting 'modified-after' time to %s" %
                self._modifiedAfterTime.toString("%FT%TZ"))
        except BaseException:
            pass

        try:
            self.evtypes = self.commandline().optionString("ev-type").split(',')
        except BaseException:
            self.evtypes = None
        try:
            self.orgType = self.commandline().optionString("org-type")
        except BaseException:
            self.orgType = "preferred"

        self.simple = self.commandline().hasOption("simple")
        self.manualOnly = self.commandline().hasOption("manual-only")
        self.automaticOnly = self.commandline().hasOption("auto-only")

        try:
            self.incAuthor = self.commandline().optionString("inc-author").split(',')
        except BaseException:
            self.incAuthor = None
        try:
            self.exclAuthor = self.commandline().optionString("excl-author").split(',')
        except BaseException:
            self.exclAuthor = None

        try:
            self.incAgencyID = self.commandline().optionString("inc-agency").split(',')
        except BaseException:
            self.incAgencyID = None

        try:
            self.exclAgencyID = self.commandline().optionString("excl-agency").split(',')
        except BaseException:
            self.exclAgencyID = None

        try:
            self.incMethodID = self.commandline().optionString("inc-method").split(',')
        except BaseException:
            self.incMethodID = None
        try:
            self.exclMethodID = self.commandline().optionString("excl-method").split(',')
        except BaseException:
            self.exclMethodID = None
        try:
            tokens = self.commandline().optionString("area").split(',')
            Area = namedtuple('Area', 'minLat minLon maxLat maxLon')
            self.area = Area(float(tokens[0]), float(tokens[1]),
                             float(tokens[2]), float(tokens[3]))
        except BaseException:
            self.area = None

        return True

    def run(self):

        if not self.simple:
            sys.stdout.write(
                "origin,event,eventType,evalMode,agencyID,author,methodID,RMS,numPhases,latitude,longitude,depth,time,creationTime,modificationTime\n")

        events = []
        for obj in self.query().getEvents(self._startTime, self._endTime):
            evt = sc_datamodel.Event.Cast(obj)
            if not evt:
                continue

            if self._modifiedAfterTime is not None:
                try:
                    if evt.creationInfo().modificationTime() < self._modifiedAfterTime:
                        continue
                except ValueError:
                    continue

            try:
                evtype = sc_datamodel.EEventTypeNames.name(evt.type())
            except ValueError:
                evtype = None
            if evtype is None:
              evtype = "" # this allows the user to match the empty event type
            if (self.evtypes is not None) and (evtype not in self.evtypes):
                continue

            events.append((evt.publicID(), evt.preferredOriginID(), evtype))

        for (evId, preferredOrgId, evtype) in events:
            orgInfo = self.findOrigin(evId, preferredOrgId)
            if orgInfo is None:
                continue
            if self.simple:
                sys.stdout.write("%s\n" % orgInfo.id)
            else:
                sys.stdout.write(
                    f"{orgInfo.id},{evId},{evtype},{orgInfo.evalMode},{orgInfo.agencyID},{orgInfo.author},{orgInfo.methodID},"
                    f"{'' if orgInfo.rms is None else orgInfo.rms},"
                    f"{'' if orgInfo.phases is None else orgInfo.phases},"
                    f"{orgInfo.latitude},{orgInfo.longitude},{orgInfo.depth},{orgInfo.time},"
                    f"{'' if orgInfo.creationTime is None else orgInfo.creationTime},"
                    f"{'' if orgInfo.modificationTime is None else orgInfo.modificationTime}\n")

        return True

    def findOrigin(self, evId, preferredOrgId):
        origins = []
        if self.orgType == "preferred":
            obj = self.query().getObject(sc_datamodel.Origin.TypeInfo(), preferredOrgId)
            org = sc_datamodel.Origin.Cast(obj)
            if org is not None:
                origins.append(org)
        else:
            for obj in self.query().getOrigins(evId):
                org = sc_datamodel.Origin.Cast(obj)
                if org is not None:
                    origins.append(org)

        OrgInfo = namedtuple(
            'OrgInfo', 'id evalMode creationTime modificationTime agencyID author methodID rms phases latitude longitude depth time sourceOrigin')
        orgInfo = None
        for currOrg in origins:

            try:
                quality = currOrg.quality()
            except ValueError:
                quality = None

            try:
                evalMode = currOrg.evaluationMode()
            except ValueError:
                evalMode = None

            try:
                creationTime = currOrg.creationInfo().creationTime()
            except ValueError:
                creationTime = None

            try:
                modificationTime = currOrg.creationInfo().modificationTime()
            except ValueError:
                modificationTime = None

            try:
                author = currOrg.creationInfo().author()
            except ValueError:
                author = None

            try:
                agencyID = currOrg.creationInfo().agencyID()
            except ValueError:
                agencyID = None

            try:
                methodID = currOrg.methodID()
            except ValueError:
                methodID = None

            try:
                latitude = currOrg.latitude().value()
                longitude = currOrg.longitude().value()
                depth = currOrg.depth().value()
            except ValueError:
                latitude = None
                longitude = None
                depth = None

            try:
                time = currOrg.time().value()
            except ValueError:
                time = None

            try:
                time = currOrg.time().value()
            except ValueError:
                time = None

            try:
              phases = quality.usedPhaseCount()
            except:
              phases = None

            try:
              rms = quality.standardError()
            except:
              rms = None

            if self.automaticOnly and evalMode != sc_datamodel.AUTOMATIC:
                continue

            if self.manualOnly and evalMode != sc_datamodel.MANUAL:
                continue

            if self.orgType == "last" and (orgInfo is not None) and \
               (creationTime < orgInfo.creationTime):
                continue

            if self.orgType == "first" and (orgInfo is not None) and \
               (creationTime > orgInfo.creationTime):
                continue

            if self.orgType == "maxPhases" and (phases is None or
               (orgInfo is not None and phases < orgInfo.phases)):
                continue

            if self.orgType == "minPhases" and (phases is None or
               (orgInfo is not None and phases > orgInfo.phases)):
                continue

            if self.orgType == "maxRMS" and (rms is None or
               (orgInfo is not None and rms < orgInfo.rms)):
                continue

            if self.orgType == "minRMS" and (rms is None or
               (orgInfo is not None and rms > orgInfo.rms)):
                continue

            if (self.incAuthor is not None) and (author is None or
               author not in self.incAuthor):
                continue

            if (self.exclAuthor is not None) and (author is not None) and \
               (author in self.exclAuthor):
                continue

            if (self.incAgencyID is not None) and (agencyID is None or
                agencyID not in self.incAgencyID):
                continue

            if (self.exclAgencyID is not None) and (agencyID is not None) and \
                    (agencyID in self.exclAgencyID):
                continue

            if (self.incMethodID is not None) and (methodID is None or
                methodID not in self.incMethodID):
                continue

            if (self.exclMethodID is not None) and (methodID is not None) and \
               (methodID in self.exclMethodID):
                continue

            if self.area is not None:

                if latitude < self.area.minLat or latitude > self.area.maxLat:
                    continue

                lonRange = self.area.maxLon - self.area.minLon
                if lonRange < 0:
                    lonRange += 360.0
                lonDelta = longitude - self.area.minLon
                if lonDelta < 0:
                    lonDelta += 360.0
                if lonDelta > lonRange:
                    continue

            evalModeStr = ""
            if evalMode is not None:
                evalModeStr = sc_datamodel.EEvaluationModeNames.name(evalMode)

            sourceOrigin = ""
            # self.query().loadComments(currOrg)
            # for i in range(currOrg.commentCount()):
            #    comment = sc_datamodel.Comment.Cast(currOrg.comment(i))
            #    if comment is None:
            #        continue
            #    if comment.id() == "relocation::sourceOrigin":
            #        sourceOrigin = comment.text()
            #        break

            orgInfo = OrgInfo(currOrg.publicID(), evalModeStr, creationTime,
                              modificationTime, agencyID, author, methodID,
                              rms, phases, latitude, longitude,  depth, time,
                              sourceOrigin)

        return orgInfo


def main():
    app = EventList(len(sys.argv), sys.argv)
    app()


if __name__ == "__main__":
    main()
