"use client"

import * as React from "react"
import { AudioWaveform, Bot, BookOpen, Command, GalleryVerticalEnd, Settings2, SquareTerminal, Frame, PieChart, Map } from 'lucide-react';

import { Separator } from "@/components/ui/separator";
import {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarGroup,
  SidebarGroupLabel,
  SidebarHeader,
  SidebarInset,
  SidebarMenu,
  SidebarMenuAction,
  SidebarMenuButton,
  SidebarMenuItem,
  SidebarMenuSub,
  SidebarMenuSubButton,
  SidebarMenuSubItem,
  SidebarProvider,
  SidebarRail,
  SidebarTrigger,
} from "@/components/ui/sidebar"
import { Dashboard } from "./dashboard"
import { UserMenu } from "./user-menu";
import { ProjectsMenu } from "./projects-menu";
import { TeamsMenu } from "./teams-menu";
import { PlatformMenu } from "./platform-menu";
// This is sample data.
const data = {
  user: {
    name: "example",
    email: "m@example.com",
    avatar: "/avatars/shadcn.jpg",
  },
  teams: [
    {
      name: "Example Corp.",
      logo: GalleryVerticalEnd,
      plan: "Enterprise",
    },
    {
      name: "Example Corp.",
      logo: AudioWaveform,
      plan: "Startup",
    },
    {
      name: "Example Corp.",
      logo: Command,
      plan: "Free",
    },
  ],
  navMain: [
    {
      title: "Analysis",
      url: "#",
      icon: SquareTerminal,
      isActive: true,
      items: [
        {
          title: "Alpha Diversity",
          url: "#",
        },
        {
          title: "Beta Diversity",
          url: "#",
        },
        {
          title: "Core Microbiome",
          url: "#",
        },
        {
          title: "Environmental Filtering",
          url: "#",
        },
        {
          title: "Community Assembly",
          url: "#",
        },
        {
          title: "Subset Regression",
          url: "#",
        },
        {
          title: "Differential Expression",
          url: "#",
        },
        {
          title: "CODA LASSO",
          url: "#",
        },
      ],
    },
    {
      title: "Documentation",
      url: "#",
      icon: BookOpen,
      items: [
        {
          title: "Introduction",
          url: "#",
        },
        {
          title: "Get Started",
          url: "#",
        },
        {
          title: "Tutorials",
          url: "#",
        },
        {
          title: "Changelog",
          url: "#",
        },
      ],
    },
    {
      title: "Settings",
      url: "#",
      icon: Settings2,
      items: [
        {
          title: "General",
          url: "#",
        },
        {
          title: "Team",
          url: "#",
        },
        {
          title: "Billing",
          url: "#",
        },
        {
          title: "Limits",
          url: "#",
        },
      ],
    },
  ],
  projects: [
    {
      name: "Example 1",
      url: "#",
      icon: Frame,
    },
    {
      name: "Example 2",
      url: "#",
      icon: PieChart,
    },
    {
      name: "Example 3",
      url: "#",
      icon: Map,
    },
  ],
}

export function SidebarLayout() {
  return (
    <SidebarProvider>
      <Sidebar collapsible="icon">
        <SidebarHeader>
          <SidebarMenu>
            <SidebarMenuItem>
              <TeamsMenu teams={data.teams} />
            </SidebarMenuItem>
          </SidebarMenu>
        </SidebarHeader>
        <SidebarContent>
          <PlatformMenu navMain={data.navMain} />
          <ProjectsMenu projects={data.projects} />
        </SidebarContent>
        <SidebarFooter>
          <SidebarMenu>
            <SidebarMenuItem>
              <UserMenu user={data.user} />
            </SidebarMenuItem>
          </SidebarMenu>
        </SidebarFooter>
        <SidebarRail />
      </Sidebar>
      <SidebarInset>
        <header className="flex h-16 shrink-0 items-center gap-2 border-b px-6">
          {/* <SidebarTrigger /> */}
          {/* <Separator orientation="vertical" className="h-6" /> */}
          <div className="flex items-center gap-5 text-sm">
            <h2 className="font-semibold">Ecology Analysis</h2>
          </div>
        </header>
        <main className="flex-1 overflow-y-auto">
          <div className="container mx-auto p-6">
            <Dashboard />
          </div>
        </main>
      </SidebarInset>
    </SidebarProvider>
  )
}
